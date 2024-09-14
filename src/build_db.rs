// build_db.rs

use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::process;
use std::sync::{Arc, Mutex};

use crate::compact_hash::CompactHashTable;

use crate::kraken2_data::IndexOptions;
use crate::kv_store::HValue;
use crate::mmscanner::{
    MinimizerScanner, CURRENT_REVCOM_VERSION, DEFAULT_SPACED_SEED_MASK, DEFAULT_TOGGLE_MASK,
};
use crate::seqreader::{BatchSequenceReader, Sequence};
use crate::taxonomy::{NCBITaxonomy, TaxId, Taxonomy};
use atty::{self, Stream};
use clap::Parser;
use libc::isatty;
use rayon::prelude::*;

const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

/// Kraken 2 Database Builder Options
#[derive(Debug, Parser)]
#[command(
    name = "Kraken 2 Database Builder",
    version = "1.0",
    about = "Builds Kraken 2 database"
)]
struct Options {
    /// Kraken 2 hash table filename
    #[arg(short = 'H')]
    hashtable_filename: String,

    /// Sequence ID to taxon map filename
    #[arg(short = 'm')]
    id_to_taxon_map_filename: String,

    /// Kraken 2 taxonomy filename
    #[arg(short = 't')]
    taxonomy_filename: String,

    /// NCBI taxonomy directory name
    #[arg(short = 'n')]
    ncbi_taxonomy_directory: String,

    /// Kraken 2 options filename
    #[arg(short = 'o')]
    options_filename: String,

    /// Set length of k-mers
    #[arg(short = 'k')]
    k: isize,

    /// Set length of minimizers
    #[arg(short = 'l')]
    l: isize,

    /// Set capacity of hash table
    #[arg(short = 'c')]
    capacity: usize,

    /// Set maximum capacity of hash table (MiniKraken)
    #[arg(short = 'M', default_value_t = 0)]
    maximum_capacity: usize,

    /// Spaced seed mask
    #[arg(short = 'S', default_value_t = DEFAULT_SPACED_SEED_MASK)]
    spaced_seed_mask: u64,

    /// Minimizer toggle mask
    #[arg(short = 'T', default_value_t = DEFAULT_TOGGLE_MASK)]
    toggle_mask: u64,

    /// Input sequences are proteins
    #[arg(short = 'X', default_value_t = false)]
    input_is_protein: bool,

    /// Number of threads
    #[arg(short = 'p', default_value_t = 1)]
    num_threads: usize,

    /// Use fast, nondeterministic building method
    #[arg(short = 'F')]
    nondeterministic_build: bool,

    /// Read block size
    #[arg(short = 'B', default_value_t = DEFAULT_BLOCK_SIZE)]
    block_size: usize,

    /// Read subblock size
    #[arg(short = 'b', default_value_t = DEFAULT_SUBBLOCK_SIZE)]
    subblock_size: usize,

    /// Bit storage requested for taxid
    #[arg(short = 'r', default_value_t = 0)]
    requested_bits_for_taxid: usize,
}

fn main() -> io::Result<()> {
    let opts = Options::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    if opts.k < 1 {
        eprintln!("k must be positive integer");
        process::exit(1);
    }
    if opts.l < 1 || opts.l > 31 {
        eprintln!("l must be positive integer and no more than 31");
        process::exit(1);
    }
    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        process::exit(1);
    }
    if opts.capacity < 1 {
        eprintln!("Capacity must be positive integer");
        process::exit(1);
    }
    if opts.block_size < opts.subblock_size {
        eprintln!("Block size cannot be less than subblock size");
        process::exit(1);
    }
    if opts.maximum_capacity > 0 && opts.maximum_capacity > opts.capacity {
        eprintln!("Maximum capacity option shouldn't specify larger capacity than normal");
        process::exit(1);
    }

    let mut id_to_taxon_map = HashMap::new();
    read_id_to_taxon_map(&mut id_to_taxon_map, &opts.id_to_taxon_map_filename)?;
    generate_taxonomy(&opts, &id_to_taxon_map)?;

    eprintln!("Taxonomy parsed and converted.");

    let mut taxonomy = Taxonomy::new(&opts.taxonomy_filename, false)?;
    taxonomy.generate_external_to_internal_id_map();
    let mut bits_needed_for_value = 1;
    while (1 << bits_needed_for_value) < taxonomy.node_count as usize {
        bits_needed_for_value += 1;
    }
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        eprintln!("More bits required for storing taxid");
        process::exit(1);
    }

    let bits_for_taxid = if bits_needed_for_value < opts.requested_bits_for_taxid {
        opts.requested_bits_for_taxid
    } else {
        bits_needed_for_value
    };

    let actual_capacity = if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            eprintln!("Maximum capacity larger than requested capacity");
            process::exit(1);
        }
        opts.maximum_capacity
    } else {
        opts.capacity
    };

    let min_clear_hash_value = if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        ((1.0 - frac) * u64::MAX as f64) as u64
    } else {
        0
    };

    let kraken_table = CompactHashTable::new(actual_capacity, 32 - bits_for_taxid, bits_for_taxid)?;
    let kraken_index = Arc::new(Mutex::new(kraken_table));

    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    let taxonomy = Arc::new(taxonomy);

    if opts.nondeterministic_build {
        process_sequences_fast(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            min_clear_hash_value,
        )?;
    } else {
        process_sequences(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            min_clear_hash_value,
        )?;
    }

    eprint!("Writing data to disk... ");
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
    let mut opts_fs = OpenOptions::new()
        .write(true)
        .create(true)
        .open(&opts.options_filename)?;
    opts_fs.write_all(&bincode::serialize(&index_opts).unwrap())?;
    opts_fs.flush()?;
    eprintln!(" complete.");

    Ok(())
}

// The rest of your code remains the same, with adjustments as per the solutions above.
// Make sure to adjust function signatures, implement Clone, and handle lifetimes as needed.
// Function to read the ID to taxon map
fn read_id_to_taxon_map(id_map: &mut HashMap<String, TaxId>, filename: &str) -> io::Result<()> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.trim().split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            let taxid: TaxId = taxid_str.parse().unwrap_or(0);
            if taxid != 0 {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
    Ok(())
}

// Function to generate the taxonomy
fn generate_taxonomy(opts: &Options, id_map: &HashMap<String, TaxId>) -> io::Result<()> {
    let nodes_filename = format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory);
    let names_filename = format!("{}/names.dmp", opts.ncbi_taxonomy_directory);

    let mut ncbi_taxonomy = NCBITaxonomy::new(&nodes_filename, &names_filename)?;
    for &taxid in id_map.values() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid);
        }
    }

    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename)?;

    Ok(())
}

// Function to process sequences deterministically
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
    sequence_reader.set_reader(Box::new(reader));
    sequence_reader.set_block_size(opts.block_size);

    let mut sequence = Sequence::new();

    while sequence_reader.next_sequence(&mut sequence) {
        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let mut taxid = 0;

        for seqid in all_sequence_ids {
            if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                if ext_taxid != 0 {
                    let internal_taxid = taxonomy.get_internal_id(ext_taxid);
                    taxid = taxonomy.lowest_common_ancestor(taxid, internal_taxid);
                }
            }
        }

        if taxid != 0 {
            // Add terminator for protein sequences if not already there
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

        if atty::is(Stream::Stderr) {
            eprint!(
                "\rProcessed {} sequences ({} {})...",
                processed_seq_ct,
                processed_ch_ct,
                if opts.input_is_protein { "aa" } else { "bp" }
            );
        }
    }

    if atty::is(Stream::Stderr) {
        eprint!("\r");
    }

    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );

    Ok(())
}

// Function to process sequences nondeterministically
fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, TaxId>,
    kraken_index: &Arc<Mutex<CompactHashTable>>,
    taxonomy: &Taxonomy,
    min_clear_hash_value: u64,
) -> io::Result<()> {
    let processed_seq_ct = Arc::new(Mutex::new(0usize));
    let processed_ch_ct = Arc::new(Mutex::new(0usize));

    let stdin = io::stdin();
    let reader = stdin.lock();

    let id_to_taxon_map = Arc::new(id_to_taxon_map.clone());
    let taxonomy = Arc::new(taxonomy);

    let mut sequence_reader = BatchSequenceReader::new(reader, opts.block_size);
    sequence_reader.set_reader(Box::new(reader));
    sequence_reader.set_block_size(opts.block_size);

    let mut sequences = Vec::new();
    let mut sequence = Sequence::new();

    while sequence_reader.next_sequence(&mut sequence) {
        sequences.push(sequence.clone());
    }

    sequences.into_par_iter().for_each(|sequence| {
        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let mut taxid = 0;

        for seqid in all_sequence_ids {
            if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                if ext_taxid != 0 {
                    let internal_taxid = taxonomy.get_internal_id(ext_taxid);
                    taxid = taxonomy.lowest_common_ancestor(taxid, internal_taxid);
                }
            }
        }

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
                CURRENT_REVCOM_VERSION,
            );

            process_sequence_fast(
                &seq,
                taxid,
                kraken_index.clone(),
                taxonomy.clone(),
                &mut scanner,
                min_clear_hash_value,
            );

            let mut seq_ct = processed_seq_ct.lock().unwrap();
            *seq_ct += 1;

            let mut ch_ct = processed_ch_ct.lock().unwrap();
            *ch_ct += seq.len();
        }
    });

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

// Function to process a single sequence deterministically
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
        // In the closure:
        let set = &minimizer_sets[zone];
        let mut set_guard = set.lock().unwrap();
        set_guard.insert(minimizer);

        // Process subblocks in parallel
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
                    CURRENT_REVCOM_VERSION,
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
                    minimizer_sets[zone].insert(minimizer);
                }
            });

        // Combine sets into a single list
        let mut minimizer_list = Vec::new();
        for set in minimizer_sets {
            let mut vec: Vec<u64> = set.into_iter().collect();
            vec.sort_unstable();
            minimizer_list.extend(vec);
        }

        // Deterministic order
        minimizer_list.sort_unstable();
        minimizer_list.dedup();

        // Set minimizer LCA in hash table
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

// Function to process a single sequence nondeterministically
fn process_sequence_fast<'a>(
    seq: &'a str,
    taxid: TaxId,
    kraken_index: Arc<Mutex<CompactHashTable>>,
    taxonomy: Arc<Taxonomy>,
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
            taxonomy.as_ref(),
        );
    }
}

// Function to extract NCBI sequence IDs from a header
fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = String::new();

    let mut in_id = true;
    let chars: Vec<char> = header.chars().collect();
    let mut i = 1; // Start after '>'
    while i < chars.len() {
        if chars[i] == '\x01' {
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = true;
        } else if in_id && chars[i].is_whitespace() {
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = false;
        } else if in_id {
            current_str.push(chars[i]);
        }
        i += 1;
    }

    if !current_str.is_empty() {
        list.push(current_str);
    }

    list
}

// Function to set the LCA for a minimizer in the hash table
fn set_minimizer_lca(
    kraken_index: &mut CompactHashTable,
    minimizer: u64,
    taxid: TaxId,
    taxonomy: &Taxonomy,
) {
    let mut existing_taxid = 0;
    let mut new_taxid = taxid;

    while !kraken_index.compare_and_set(
        minimizer,
        (new_taxid as HValue).into(),
        &mut existing_taxid,
    ) {
        new_taxid = taxonomy.lowest_common_ancestor(new_taxid, existing_taxid as TaxId);
    }
}

// Placeholder for murmur_hash3 function
fn murmur_hash3(key: u64) -> u64 {
    // Implement MurmurHash3 here
    key // Placeholder implementation
}
