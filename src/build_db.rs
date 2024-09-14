// build_db.rs

use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::process;
use std::sync::{Arc, Mutex};

use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::compact_hash::{CompactHashTable, IndexOptions};
use crate::kraken2_data::{CURRENT_REVCOM_VERSION, DEFAULT_SPACED_SEED_MASK, DEFAULT_TOGGLE_MASK};
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{BatchSequenceReader, Sequence};
use crate::taxonomy::{NCBITaxonomy, TaxId, Taxonomy};

const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

/// Kraken 2 Database Builder
#[derive(Debug, Parser)]
#[clap(version = "1.0", about = "Kraken 2 Database Builder")]
struct Options {
    /// Kraken 2 hash table filename
    #[clap(short = 'H', long)]
    hashtable_filename: String,

    /// Sequence ID to taxon map filename
    #[clap(short = 'm', long)]
    id_to_taxon_map_filename: String,

    /// Kraken 2 taxonomy filename
    #[clap(short = 't', long)]
    taxonomy_filename: String,

    /// NCBI taxonomy directory name
    #[clap(short = 'n', long)]
    ncbi_taxonomy_directory: String,

    /// Kraken 2 options filename
    #[clap(short = 'o', long)]
    options_filename: String,

    /// Set length of k-mers
    #[clap(short = 'k', long)]
    k: usize,

    /// Set length of minimizers
    #[clap(short = 'l', long)]
    l: usize,

    /// Set capacity of hash table
    #[clap(short = 'c', long)]
    capacity: usize,

    /// Set maximum capacity of hash table (MiniKraken)
    #[clap(short = 'M', long)]
    maximum_capacity: Option<usize>,

    /// Spaced seed mask
    #[clap(short = 'S', long)]
    spaced_seed_mask: Option<String>,

    /// Minimizer toggle mask
    #[clap(short = 'T', long)]
    toggle_mask: Option<String>,

    /// Input sequences are proteins
    #[clap(short = 'X', long)]
    input_is_protein: bool,

    /// Number of threads
    #[clap(short = 'p', long, default_value = "1")]
    num_threads: usize,

    /// Use fast, nondeterministic building method
    #[clap(short = 'F', long)]
    nondeterministic_build: bool,

    /// Read block size
    #[clap(short = 'B', long, default_value = "10485760")]
    block_size: usize,

    /// Read subblock size
    #[clap(short = 'b', long, default_value = "1024")]
    subblock_size: usize,

    /// Bit storage requested for taxid
    #[clap(short = 'r', long, default_value = "0")]
    requested_bits_for_taxid: usize,
}

fn main() -> io::Result<()> {
    let opts = Options::parse();

    // Set the number of threads for Rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    // Default values for spaced seed mask and toggle mask
    let spaced_seed_mask = if let Some(s) = opts.spaced_seed_mask {
        u64::from_str_radix(&s, 2).unwrap_or(DEFAULT_SPACED_SEED_MASK)
    } else {
        DEFAULT_SPACED_SEED_MASK
    };

    let toggle_mask = if let Some(s) = opts.toggle_mask {
        u64::from_str_radix(&s, 2).unwrap_or(DEFAULT_TOGGLE_MASK)
    } else {
        DEFAULT_TOGGLE_MASK
    };

    // Validate options
    if opts.k < opts.l {
        eprintln!("Error: k cannot be less than l");
        process::exit(1);
    }

    if opts.block_size < opts.subblock_size {
        eprintln!("Error: Block size cannot be less than subblock size");
        process::exit(1);
    }

    if let Some(max_capacity) = opts.maximum_capacity {
        if max_capacity > opts.capacity {
            eprintln!("Error: Maximum capacity cannot be larger than capacity");
            process::exit(1);
        }
    }

    // Read ID to taxon map
    let id_to_taxon_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;

    // Generate taxonomy
    generate_taxonomy(&opts, &id_to_taxon_map)?;

    eprintln!("Taxonomy parsed and converted.");

    // Load taxonomy
    let mut taxonomy = Taxonomy::new(&opts.taxonomy_filename)?;
    taxonomy.generate_external_to_internal_id_map();

    let mut bits_needed_for_value = 1;
    while (1 << bits_needed_for_value) < taxonomy.node_count() {
        bits_needed_for_value += 1;
    }

    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        eprintln!("Error: More bits required for storing taxid");
        process::exit(1);
    }

    let bits_for_taxid = std::cmp::max(bits_needed_for_value, opts.requested_bits_for_taxid);

    let actual_capacity = if let Some(max_capacity) = opts.maximum_capacity {
        max_capacity
    } else {
        opts.capacity
    };

    let min_clear_hash_value = if let Some(max_capacity) = opts.maximum_capacity {
        let frac = max_capacity as f64 / opts.capacity as f64;
        ((1.0 - frac) * u64::MAX as f64) as u64
    } else {
        0
    };

    let kraken_index = Arc::new(Mutex::new(CompactHashTable::new(
        actual_capacity,
        32 - bits_for_taxid,
        bits_for_taxid,
    )?));

    eprintln!(
        "CompactHashTable created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    if !opts.nondeterministic_build {
        process_sequences(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            spaced_seed_mask,
            toggle_mask,
            min_clear_hash_value,
        )?;
    } else {
        process_sequences_fast(
            &opts,
            &id_to_taxon_map,
            &kraken_index,
            &taxonomy,
            spaced_seed_mask,
            toggle_mask,
            min_clear_hash_value,
        )?;
    }

    eprintln!("Writing data to disk... ");
    kraken_index
        .lock()
        .unwrap()
        .write_table(&opts.hashtable_filename)?;

    // Write index options
    let index_opts = IndexOptions {
        k: opts.k as u32,
        l: opts.l as u32,
        spaced_seed_mask,
        toggle_mask,
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

    eprintln!("Complete.");

    Ok(())
}

// Function to read the ID to taxon map
fn read_id_to_taxon_map(filename: &str) -> io::Result<HashMap<String, TaxId>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();

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

    Ok(id_map)
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
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
) -> io::Result<()> {
    let mut processed_seq_ct = 0;
    let mut processed_ch_ct = 0;

    let stdin = io::stdin();
    let reader = stdin.lock();

    let mut sequence_reader = BatchSequenceReader::new(reader, opts.block_size);

    while let Some(sequence) = sequence_reader.next_sequence()? {
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
                spaced_seed_mask,
                toggle_mask,
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
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
) -> io::Result<()> {
    let processed_seq_ct = Arc::new(Mutex::new(0usize));
    let processed_ch_ct = Arc::new(Mutex::new(0usize));

    let stdin = io::stdin();
    let reader = stdin.lock();

    let id_to_taxon_map = Arc::new(id_to_taxon_map.clone());
    let taxonomy = Arc::new(taxonomy.clone());

    let sequence_reader = BatchSequenceReader::new(reader, opts.block_size);

    let sequences = sequence_reader.read_all_sequences()?;

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

            process_sequence_fast(
                &seq,
                taxid,
                kraken_index.clone(),
                taxonomy.clone(),
                opts.k,
                opts.l,
                spaced_seed_mask,
                toggle_mask,
                min_clear_hash_value,
                opts.input_is_protein,
            );

            let mut seq_ct = processed_seq_ct.lock().unwrap();
            *seq_ct += 1;
            drop(seq_ct);

            let mut ch_ct = processed_ch_ct.lock().unwrap();
            *ch_ct += seq.len();
            drop(ch_ct);
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
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
) {
    let set_ct = 256;
    let mut locks = vec![Mutex::new(()); set_ct];

    for j in (0..seq.len()).step_by(opts.block_size) {
        let block_start = j;
        let block_finish = std::cmp::min(j + opts.block_size + opts.k - 1, seq.len());

        let mut minimizer_sets: Vec<HashSet<u64>> = vec![HashSet::new(); set_ct];

        // Process subblocks in parallel
        (block_start..block_finish)
            .into_par_iter()
            .step_by(opts.subblock_size)
            .for_each(|i| {
                let subblock_finish =
                    std::cmp::min(i + opts.subblock_size + opts.k - 1, block_finish);

                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    spaced_seed_mask,
                    !opts.input_is_protein,
                    toggle_mask,
                );
                scanner.load_sequence(&seq[i..subblock_finish]);

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
fn process_sequence_fast(
    seq: &str,
    taxid: TaxId,
    kraken_index: Arc<Mutex<CompactHashTable>>,
    taxonomy: Arc<Taxonomy>,
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
    input_is_protein: bool,
) {
    let mut scanner = MinimizerScanner::new(k, l, spaced_seed_mask, !input_is_protein, toggle_mask);
    scanner.load_sequence(seq);

    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        let hc = murmur_hash3(minimizer);
        if min_clear_hash_value != 0 && hc < min_clear_hash_value {
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

    while !kraken_index.compare_and_set(minimizer, new_taxid as u64, &mut existing_taxid) {
        new_taxid = taxonomy.lowest_common_ancestor(new_taxid, existing_taxid as TaxId);
    }
}

// Placeholder for murmur_hash3 function
fn murmur_hash3(key: u64) -> u64 {
    // Implement MurmurHash3 here
    // This is a placeholder implementation
    key
}
