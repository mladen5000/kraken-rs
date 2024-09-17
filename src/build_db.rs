use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;

use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::IndexOptions;
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{BatchSequenceReader, Sequence};
use crate::taxonomy::{NCBITaxonomy, TaxId, Taxonomy};

const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Options {
    #[arg(short = 'H')]
    hashtable_filename: String,

    #[arg(short = 'm')]
    id_to_taxon_map_filename: String,

    #[arg(short = 'n')]
    ncbi_taxonomy_directory: String,

    #[arg(short = 'o')]
    options_filename: String,

    #[arg(short = 't')]
    taxonomy_filename: String,

    #[arg(short = 'B', default_value_t = DEFAULT_BLOCK_SIZE)]
    block_size: usize,

    #[arg(short = 'b', default_value_t = DEFAULT_SUBBLOCK_SIZE)]
    subblock_size: usize,

    #[arg(short = 'r', default_value_t = 0)]
    requested_bits_for_taxid: usize,

    #[arg(short = 'p', default_value_t = 1)]
    num_threads: usize,

    #[arg(short = 'X')]
    input_is_protein: bool,

    #[arg(short = 'k')]
    k: isize,

    #[arg(short = 'l')]
    l: isize,

    #[arg(short = 'c')]
    capacity: usize,

    #[arg(short = 'M', default_value_t = 0)]
    maximum_capacity: usize,

    #[arg(short = 'S', default_value_t = DEFAULT_SPACED_SEED_MASK)]
    spaced_seed_mask: u64,

    #[arg(short = 'T', default_value_t = DEFAULT_TOGGLE_MASK)]
    toggle_mask: u64,

    #[arg(short = 'F')]
    deterministic_build: bool,
}

fn main() -> io::Result<()> {
    let opts = Options::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    let mut id_to_taxon_map = HashMap::new();
    read_id_to_taxon_map(&mut id_to_taxon_map, &opts.id_to_taxon_map_filename)?;
    generate_taxonomy(&opts, &id_to_taxon_map)?;

    eprintln!("Taxonomy parsed and converted.");

    let mut taxonomy = Taxonomy::new(&opts.taxonomy_filename, false)?;
    taxonomy.generate_external_to_internal_id_map();

    let bits_needed_for_value = (taxonomy.node_count as f64).log2().ceil() as usize;
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        eprintln!("More bits required for storing taxid");
        std::process::exit(1);
    }

    let bits_for_taxid = std::cmp::max(bits_needed_for_value, opts.requested_bits_for_taxid);

    let (actual_capacity, min_clear_hash_value) = if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            eprintln!("Maximum capacity larger than requested capacity");
            std::process::exit(1);
        }
        let min_clear_hash_value = ((1.0 - frac) * u64::MAX as f64) as u64;
        (opts.maximum_capacity, min_clear_hash_value)
    } else {
        (opts.capacity, 0)
    };

    let kraken_index = Arc::new(Mutex::new(CompactHashTable::new(
        actual_capacity,
        32 - bits_for_taxid,
        bits_for_taxid,
    )?));

    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    let start_time = Instant::now();

    if opts.deterministic_build {
        process_sequences(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            min_clear_hash_value,
        )?;
    } else {
        process_sequences_fast(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            min_clear_hash_value,
        )?;
    }

    eprintln!("Writing data to disk... ");
    kraken_index
        .lock()
        .unwrap()
        .write_table(&opts.hashtable_filename)?;

    let index_opts = IndexOptions {
        k: opts.k as usize,
        l: opts.l as usize,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        minimum_acceptable_hash_value: min_clear_hash_value,
        revcom_version: CURRENT_REVCOM_VERSION,
        db_version: 0,
        db_type: 0,
    };

    let mut opts_file = File::create(&opts.options_filename)?;
    opts_file.write_all(&bincode::serialize(&index_opts).unwrap())?;

    eprintln!(" complete.");
    eprintln!("Total build time: {:?}", start_time.elapsed());

    Ok(())
}

fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, TaxId>,
    kraken_index: &Arc<Mutex<CompactHashTable>>,
    taxonomy: &Taxonomy,
    min_clear_hash_value: u64,
) -> io::Result<()> {
    let mut processed_seq_ct = 0;
    let mut processed_ch_ct = 0;

    let stdin = io::stdin();
    let reader = stdin.lock();
    let mut sequence_reader = BatchSequenceReader::new(reader, opts.block_size);

    let mut sequence = Sequence::default();
    while sequence_reader.next_sequence(&mut sequence) {
        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let taxid = all_sequence_ids
            .iter()
            .filter_map(|seqid| id_to_taxon_map.get(seqid))
            .fold(0, |acc, &ext_taxid| {
                taxonomy.lowest_common_ancestor(acc, taxonomy.get_internal_id(ext_taxid))
            });

        if taxid != 0 {
            let seq = if opts.input_is_protein && !sequence.seq.ends_with('*') {
                format!("{}*", sequence.seq)
            } else {
                sequence.seq.clone()
            };

            process_sequence(
                opts,
                &seq,
                taxid,
                kraken_index,
                taxonomy,
                min_clear_hash_value,
            );

            processed_seq_ct += 1;
            processed_ch_ct += seq.len();
        }

        if atty::is(atty::Stream::Stderr) {
            eprint!(
                "\rProcessed {} sequences ({} {})...",
                processed_seq_ct,
                processed_ch_ct,
                if opts.input_is_protein { "aa" } else { "bp" }
            );
        }
    }

    if atty::is(atty::Stream::Stderr) {
        eprintln!();
    }
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );

    Ok(())
}

fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, TaxId>,
    kraken_index: &Arc<Mutex<CompactHashTable>>,
    taxonomy: &Taxonomy,
    min_clear_hash_value: u64,
) -> io::Result<()> {
    let processed_seq_ct = Arc::new(Mutex::new(0));
    let processed_ch_ct = Arc::new(Mutex::new(0));

    let stdin = io::stdin();
    let reader = stdin.lock();
    let mut sequence_reader = BatchSequenceReader::new(reader, opts.block_size);

    let mut sequence = Sequence::default();
    while sequence_reader.next_sequence(&mut sequence) {
        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let taxid = all_sequence_ids
            .iter()
            .filter_map(|seqid| id_to_taxon_map.get(seqid))
            .fold(0, |acc, &ext_taxid| {
                taxonomy.lowest_common_ancestor(acc, taxonomy.get_internal_id(ext_taxid))
            });

        if taxid != 0 {
            let seq = if opts.input_is_protein && !sequence.seq.ends_with('*') {
                format!("{}*", sequence.seq)
            } else {
                sequence.seq.clone()
            };

            let mut scanner = MinimizerScanner::new(
                opts.k,
                opts.l,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
                CURRENT_REVCOM_VERSION as i32,
            );

            process_sequence_fast(
                &seq,
                taxid,
                kraken_index,
                taxonomy,
                &mut scanner,
                min_clear_hash_value,
            );

            let mut seq_ct = processed_seq_ct.lock().unwrap();
            *seq_ct += 1;

            let mut ch_ct = processed_ch_ct.lock().unwrap();
            *ch_ct += seq.len();
        }

        if atty::is(atty::Stream::Stderr) {
            let seq_ct = *processed_seq_ct.lock().unwrap();
            let ch_ct = *processed_ch_ct.lock().unwrap();
            eprint!(
                "\rProcessed {} sequences ({} {})...",
                seq_ct,
                ch_ct,
                if opts.input_is_protein { "aa" } else { "bp" }
            );
        }
    }

    if atty::is(atty::Stream::Stderr) {
        eprintln!();
    }

    let seq_ct = *processed_seq_ct.lock().unwrap();
    let ch_ct = *processed_ch_ct.lock().unwrap();
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        seq_ct,
        ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );

    Ok(())
}

fn process_sequence(
    opts: &Options,
    seq: &str,
    taxid: TaxId,
    kraken_index: &Arc<Mutex<CompactHashTable>>,
    taxonomy: &Taxonomy,
    min_clear_hash_value: u64,
) {
    let set_ct = 256;
    let locks: Vec<Mutex<()>> = (0..set_ct).map(|_| Mutex::new(())).collect();

    for j in (0..seq.len()).step_by(opts.block_size) {
        let block_start = j;
        let block_finish = std::cmp::min(j + opts.block_size + opts.k as usize - 1, seq.len());

        let minimizer_sets: Vec<Mutex<HashSet<u64>>> =
            (0..set_ct).map(|_| Mutex::new(HashSet::new())).collect();

        (block_start..block_finish)
            .into_par_iter()
            .step_by(opts.subblock_size)
            .for_each(|i| {
                let subblock_finish =
                    std::cmp::min(i + opts.subblock_size + opts.k as usize - 1, block_finish);
                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    opts.spaced_seed_mask,
                    !opts.input_is_protein,
                    opts.toggle_mask,
                    CURRENT_REVCOM_VERSION as i32,
                );

                scanner.load_sequence(seq, i, subblock_finish);

                while let Some(minimizer) = scanner.next_minimizer() {
                    if scanner.is_ambiguous() {
                        continue;
                    }
                    let hc = murmur_hash3(minimizer);
                    if min_clear_hash_value != 0 && hc < min_clear_hash_value {
                        continue;
                    }
                    let zone = (hc % set_ct as u64) as usize;
                    let _lock = locks[zone].lock().unwrap();
                    minimizer_sets[zone].lock().unwrap().insert(minimizer);
                }
            });

        let mut minimizer_list = Vec::new();
        for set in &minimizer_sets {
            minimizer_list.extend(set.lock().unwrap().iter().cloned());
        }

        minimizer_list.sort_unstable();
        minimizer_list.dedup();

        for minimizer in minimizer_list {
            set_minimizer_lca(
                &mut kraken_index.lock().unwrap(),
                minimizer,
                taxid,
                taxonomy,
            );
        }
    }
}

fn process_sequence_fast<'a>(
    seq: &'a str,
    taxid: TaxId,
    kraken_index: &Arc<Mutex<CompactHashTable>>,
    taxonomy: &Taxonomy,
    scanner: &mut MinimizerScanner<'a>,
    min_clear_hash_value: u64,
) {
    scanner.load_sequence(seq, 0, seq.len());

    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0 && murmur_hash3(minimizer) < min_clear_hash_value {
            continue;
        }
        set_minimizer_lca(
            &mut kraken_index.lock().unwrap(),
            minimizer,
            taxid,
            taxonomy,
        );
    }
}

fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    for (i, c) in header.chars().enumerate() {
        if i == 0 {
            continue; // Skip the '>' character
        }

        if c == '\x01' {
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = true;
        } else if in_id && c.is_whitespace() {
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = false;
        } else if in_id {
            current_str.push(c);
        }
    }

    if !current_str.is_empty() {
        list.push(current_str);
    }

    list
}

fn set_minimizer_lca(
    kraken_index: &mut CompactHashTable,
    minimizer: u64,
    taxid: TaxId,
    taxonomy: &Taxonomy,
) {
    let mut existing_taxid = 0;
    let mut new_taxid = taxid;

    while !kraken_index.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
        new_taxid = taxonomy.lowest_common_ancestor(new_taxid, existing_taxid as TaxId);
    }
}

fn read_id_to_taxon_map(id_map: &mut HashMap<String, TaxId>, filename: &str) -> io::Result<()> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            let taxid: TaxId = taxid_str.parse().unwrap_or(0);
            if taxid != 0 {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
    Ok(())
}

fn generate_taxonomy(opts: &Options, id_map: &HashMap<String, TaxId>) -> io::Result<()> {
    let nodes_filename = PathBuf::from(&opts.ncbi_taxonomy_directory).join("nodes.dmp");
    let names_filename = PathBuf::from(&opts.ncbi_taxonomy_directory).join("names.dmp");

    let mut ncbi_taxonomy = NCBITaxonomy::new(
        nodes_filename.to_str().unwrap(),
        names_filename.to_str().unwrap(),
    )?;
    for &taxid in id_map.values() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid);
        }
    }

    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename)?;

    Ok(())
}

fn murmur_hash3(key: u64) -> u64 {
    // Implement MurmurHash3 here
    // This is a placeholder implementation
    key
}

// Constants
const DEFAULT_SPACED_SEED_MASK: u64 = 0x3FFFFFFFFC000000;
const DEFAULT_TOGGLE_MASK: u64 = 0xAA55AA55AA55AA55;
const CURRENT_REVCOM_VERSION: u32 = 1;

// Additional helper functions

fn get_rank_code(rank: &str) -> Option<char> {
    match rank {
        "superkingdom" => Some('D'),
        "kingdom" => Some('K'),
        "phylum" => Some('P'),
        "class" => Some('C'),
        "order" => Some('O'),
        "family" => Some('F'),
        "genus" => Some('G'),
        "species" => Some('S'),
        _ => None,
    }
}

// Error handling function
fn handle_error<T, E: std::fmt::Display>(result: Result<T, E>, message: &str) -> T {
    match result {
        Ok(value) => value,
        Err(error) => {
            eprintln!("{}: {}", message, error);
            std::process::exit(1);
        }
    }
}

// Function to validate command-line arguments
fn validate_arguments(opts: &Options) {
    if opts.k < 1 {
        eprintln!("k must be a positive integer");
        std::process::exit(1);
    }
    if opts.l < 1 || opts.l > 31 {
        eprintln!("l must be a positive integer and no more than 31");
        std::process::exit(1);
    }
    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        std::process::exit(1);
    }
    if opts.capacity < 1 {
        eprintln!("Capacity must be a positive integer");
        std::process::exit(1);
    }
    if opts.block_size < opts.subblock_size {
        eprintln!("Block size cannot be less than subblock size");
        std::process::exit(1);
    }
    if opts.maximum_capacity > 0 && opts.maximum_capacity > opts.capacity {
        eprintln!("Maximum capacity option shouldn't specify larger capacity than normal");
        std::process::exit(1);
    }
}

// Main function with error handling
fn run() -> io::Result<()> {
    let opts = Options::parse();
    validate_arguments(&opts);

    // ... (rest of the main function logic)

    Ok(())
}
