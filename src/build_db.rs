//
// Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
//
// This file is part of the Kraken 2 taxonomic sequence classification system.
//

//! This Rust module is the translated version of `build_db.cc` from Kraken 2.
//! It handles command-line parsing, reading sequences, building a compact hash table,
//! and finalizing the Kraken 2 database building process.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::process::exit;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::{bail, Context, Result};
use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;

use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::IndexOptions;
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{BatchSequenceReader, Sequence};
use crate::taxonomy::{NCBITaxonomyOps, Taxonomy};
use crate::utilities::murmur_hash3_64;

pub const DEFAULT_SPACED_SEED_MASK: u64 = 0xAAAA;
pub const DEFAULT_TOGGLE_MASK: u64 = 0x5555;
pub const CURRENT_REVCOM_VERSION: u32 = 1;

// The default block sizes from the C++ code
const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

// We'll define Options as a Rust struct to match the original one.
#[derive(Debug)]
pub struct Options {
    pub id_to_taxon_map_filename: String,
    pub ncbi_taxonomy_directory: String,
    pub hashtable_filename: String,
    pub options_filename: String,
    pub taxonomy_filename: String,
    pub block_size: usize,
    pub subblock_size: usize,
    pub requested_bits_for_taxid: usize,
    pub num_threads: usize,
    pub input_is_protein: bool,
    pub k: usize,
    pub l: usize,
    pub capacity: usize,
    pub maximum_capacity: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub min_clear_hash_value: u64,
    pub deterministic_build: bool,
}

// We define a usage() function analogous to the C++ usage().
fn usage(exit_code: i32) -> ! {
    eprintln!(
        "Usage: build_db <options>\n
Options (*mandatory):
* -H FILENAME   Kraken 2 hash table filename
* -m FILENAME   Sequence ID to taxon map filename
* -t FILENAME   Kraken 2 taxonomy filename
* -n DIR        NCBI taxonomy directory name
* -o FILENAME   Kraken 2 options filename
* -k INT        Set length of k-mers
* -l INT        Set length of minimizers
* -c INT        Set capacity of hash table
  -M INT        Set maximum capacity of hash table (MiniKraken)
  -S BITSTRING  Spaced seed mask
  -T BITSTRING  Minimizer toggle mask
  -X            Input seqs. are proteins
  -p INT        Number of threads
  -F            Use fast, nondeterministic building method
  -B INT        Read block size
  -b INT        Read subblock size
  -r INT        Bit storage requested for taxid
"
    );
    exit(exit_code);
}

// Convert lines like in the C++ ParseCommandLine, but we’ll use Clap or a hand-rolled parser
fn parse_command_line() -> Result<Options> {
    let matches = Command::new("build_db")
        .about("Build the Kraken 2 database in Rust")
        .arg(
            Arg::new("H")
                .long("hashtable")
                .short('H')
                .value_name("FILENAME")
                .required(false)
                .help("Kraken 2 hash table filename")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("m")
                .long("map")
                .short('m')
                .value_name("FILENAME")
                .required(false)
                .help("Sequence ID to taxon map filename")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("t")
                .long("taxonomy")
                .short('t')
                .value_name("FILENAME")
                .required(false)
                .help("Kraken 2 taxonomy filename")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("n")
                .long("ncbi")
                .short('n')
                .value_name("DIR")
                .required(false)
                .help("NCBI taxonomy directory name")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("o")
                .long("options")
                .short('o')
                .value_name("FILENAME")
                .required(false)
                .help("Kraken 2 options filename")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("k")
                .long("klen")
                .short('k')
                .value_name("INT")
                .required(false)
                .help("Set length of k-mers")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("l")
                .long("mlen")
                .short('l')
                .value_name("INT")
                .required(false)
                .help("Set length of minimizers")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("c")
                .long("capacity")
                .short('c')
                .value_name("INT")
                .required(false)
                .help("Set capacity of hash table")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("M")
                .long("maxcap")
                .short('M')
                .value_name("INT")
                .required(false)
                .help("Set maximum capacity of hash table (MiniKraken)")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("S")
                .long("spaced")
                .short('S')
                .value_name("BITSTRING")
                .required(false)
                .help("Spaced seed mask")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("T")
                .long("toggle")
                .short('T')
                .value_name("BITSTRING")
                .required(false)
                .help("Minimizer toggle mask")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("X")
                .long("protein")
                .short('X')
                .required(false)
                .help("Input seqs. are proteins")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("p")
                .long("threads")
                .short('p')
                .value_name("INT")
                .required(false)
                .help("Number of threads")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("F")
                .long("fast")
                .short('F')
                .required(false)
                .help("Use fast, nondeterministic building method")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("B")
                .long("block")
                .short('B')
                .value_name("INT")
                .required(false)
                .help("Read block size")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("b")
                .long("subblock")
                .short('b')
                .value_name("INT")
                .required(false)
                .help("Read subblock size")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("r")
                .long("requested-bits")
                .short('r')
                .value_name("INT")
                .required(false)
                .help("Bit storage requested for taxid")
                .action(ArgAction::Set),
        )
        .get_matches();

    // Extract the matches or use defaults
    let hashtable_filename = matches
        .get_one::<String>("H")
        .unwrap_or(&"".to_string())
        .clone();
    let id_to_taxon_map_filename = matches
        .get_one::<String>("m")
        .unwrap_or(&"".to_string())
        .clone();
    let taxonomy_filename = matches
        .get_one::<String>("t")
        .unwrap_or(&"".to_string())
        .clone();
    let ncbi_taxonomy_directory = matches
        .get_one::<String>("n")
        .unwrap_or(&"".to_string())
        .clone();
    let options_filename = matches
        .get_one::<String>("o")
        .unwrap_or(&"".to_string())
        .clone();
    let k = matches
        .get_one::<String>("k")
        .map_or(0, |s| s.parse::<usize>().unwrap_or(0));
    let l = matches
        .get_one::<String>("l")
        .map_or(0, |s| s.parse::<usize>().unwrap_or(0));
    let capacity = matches
        .get_one::<String>("c")
        .map_or(0, |s| s.parse::<usize>().unwrap_or(0));
    let maximum_capacity = matches
        .get_one::<String>("M")
        .map_or(0, |s| s.parse::<usize>().unwrap_or(0));
    let spaced_seed_mask_str = matches
        .get_one::<String>("S")
        .map_or(None, |s| Some(s.clone()));
    let toggle_mask_str = matches
        .get_one::<String>("T")
        .map_or(None, |s| Some(s.clone()));
    let input_is_protein = matches.get_flag("X");
    let num_threads = matches
        .get_one::<String>("p")
        .map_or(1, |s| s.parse::<usize>().unwrap_or(1));
    let deterministic_build = !matches.get_flag("F");
    let block_size = matches
        .get_one::<String>("B")
        .map_or(DEFAULT_BLOCK_SIZE, |s| {
            s.parse::<usize>().unwrap_or(DEFAULT_BLOCK_SIZE)
        });
    let subblock_size = matches
        .get_one::<String>("b")
        .map_or(DEFAULT_SUBBLOCK_SIZE, |s| {
            s.parse::<usize>().unwrap_or(DEFAULT_SUBBLOCK_SIZE)
        });
    let requested_bits_for_taxid = matches
        .get_one::<String>("r")
        .map_or(0, |s| s.parse::<usize>().unwrap_or(0));

    // Convert spaced_seed_mask and toggle_mask from bitstring to u64
    let mut spaced_seed_mask = if let Some(s) = spaced_seed_mask_str {
        u64::from_str_radix(&s, 2).unwrap_or(DEFAULT_SPACED_SEED_MASK)
    } else {
        DEFAULT_SPACED_SEED_MASK
    };

    let toggle_mask = if let Some(s) = toggle_mask_str {
        u64::from_str_radix(&s, 2).unwrap_or(DEFAULT_TOGGLE_MASK)
    } else {
        DEFAULT_TOGGLE_MASK
    };

    if hashtable_filename.is_empty()
        || id_to_taxon_map_filename.is_empty()
        || ncbi_taxonomy_directory.is_empty()
        || options_filename.is_empty()
        || taxonomy_filename.is_empty()
    {
        eprintln!("missing mandatory filename parameter");
        usage(1);
    }

    if k == 0 || l == 0 || capacity == 0 {
        eprintln!("missing mandatory integer parameter");
        usage(1);
    }

    if k < l as usize {
        eprintln!("k cannot be less than l");
        usage(1);
    }
    if block_size < subblock_size {
        eprintln!("block size cannot be less than subblock size");
        usage(1);
    }
    if maximum_capacity > capacity {
        eprintln!("maximum capacity option shouldn't specify larger capacity than normal");
        usage(1);
    }

    if spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        let bits_per_char = if input_is_protein {
            crate::kraken2_data::BITS_PER_CHAR_PRO
        } else {
            crate::kraken2_data::BITS_PER_CHAR_DNA
        };
        expand_spaced_seed_mask(&mut spaced_seed_mask, bits_per_char);
    }

    if num_threads > rayon::current_num_threads() {
        bail!(
            "OMP only wants you to use {} threads",
            rayon::current_num_threads()
        );
    }

    Ok(Options {
        id_to_taxon_map_filename,
        ncbi_taxonomy_directory,
        hashtable_filename,
        options_filename,
        taxonomy_filename,
        block_size,
        subblock_size,
        requested_bits_for_taxid,
        num_threads,
        input_is_protein,
        k,
        l,
        capacity,
        maximum_capacity,
        spaced_seed_mask,
        toggle_mask,
        min_clear_hash_value: 0,
        deterministic_build,
    })
}

fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    for c in header.chars().skip(1) {
        // Skip '>' character
        if c == '\x01' {
            if !current_str.is_empty() {
                list.push(current_str.clone());
            }
            current_str.clear();
            in_id = true;
        } else if in_id && c.is_whitespace() {
            if !current_str.is_empty() {
                list.push(current_str.clone());
            }
            current_str.clear();
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

fn read_id_to_taxon_map(filename: &str) -> Result<HashMap<String, u64>> {
    let file =
        File::open(filename).with_context(|| format!("unable to read from '{}'", filename))?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            let taxid = taxid_str.parse::<u64>().unwrap_or(0);
            if taxid != 0 {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
    Ok(id_map)
}

fn generate_taxonomy(opts: &mut Options, id_map: &HashMap<String, u64>) -> Result<()> {
    // Build NCBI taxonomy from nodes.dmp and names.dmp, then convert to Kraken
    let ncbi_taxonomy = Taxonomy::from_ncbi_dmp(
        &format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory),
        &format!("{}/names.dmp", opts.ncbi_taxonomy_directory),
    )?;

    for taxid in id_map.values() {
        if *taxid != 0 {
            NCBITaxonomyOps::mark_node(&ncbi_taxonomy, *taxid);
        }
    }
    NCBITaxonomyOps::convert_to_kraken_taxonomy(&ncbi_taxonomy, &opts.taxonomy_filename)?;
    Ok(())
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
struct MinimizerTaxidPair {
    minimizer: u64,
    taxid: u64,
}

// A quick but nondeterministic build
fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, u64>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &mut Taxonomy,
) -> Result<()> {
    let processed_seq_ct = Arc::new(AtomicUsize::new(0));
    let processed_ch_ct = Arc::new(AtomicUsize::new(0));
    let taxonomy = Arc::new(Mutex::new(taxonomy));
    let kraken_index = Arc::new(Mutex::new(kraken_index));

    taxonomy
        .lock()
        .unwrap()
        .generate_external_to_internal_id_map();

    let stdin = io::stdin();
    let mut reader = BatchSequenceReader::new(stdin.lock());
    let mut sequences = Vec::new();

    while reader.load_block(opts.block_size)? {
        let mut sequence = Sequence::default();
        while reader.next_sequence(&mut sequence)? {
            sequences.push(sequence.clone());
        }
    }

    rayon::scope(|s| {
        for _ in 0..opts.num_threads {
            let sequences = sequences.clone();
            let processed_seq_ct = Arc::clone(&processed_seq_ct);
            let processed_ch_ct = Arc::clone(&processed_ch_ct);
            let taxonomy = Arc::clone(&taxonomy);
            let kraken_index = Arc::clone(&kraken_index);
            let opts = opts.clone();
            let id_to_taxon_map = id_to_taxon_map.clone();

            s.spawn(move |_| {
                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    opts.spaced_seed_mask,
                    !opts.input_is_protein,
                    opts.toggle_mask,
                    true,
                );

                for sequence in &sequences {
                    let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
                    let mut taxid = 0;

                    {
                        let taxonomy = taxonomy.lock().unwrap();
                        for seqid in &all_sequence_ids {
                            if let Some(&ext_taxid) = id_to_taxon_map.get(seqid) {
                                if ext_taxid != 0 {
                                    taxid = taxonomy.lowest_common_ancestor(
                                        taxid,
                                        taxonomy.get_internal_id(ext_taxid),
                                    );
                                }
                            }
                        }
                    }

                    if taxid != 0 {
                        let mut seq = sequence.seq.clone();
                        if opts.input_is_protein && !seq.ends_with('*') {
                            seq.push('*');
                        }

                        scanner.load_sequence(&seq, 0, seq.len());

                        while let Some(minimizer) = scanner.next_minimizer() {
                            if !scanner.is_ambiguous() {
                                if opts.min_clear_hash_value == 0
                                    || murmur_hash3_64(&minimizer.to_le_bytes())
                                        >= opts.min_clear_hash_value
                                {
                                    let mut hash = kraken_index.lock().unwrap();
                                    let taxonomy = taxonomy.lock().unwrap();
                                    hash.set_minimizer_lca(minimizer, taxid, &taxonomy);
                                }
                            }
                        }

                        processed_seq_ct.fetch_add(1, Ordering::Relaxed);
                        processed_ch_ct.fetch_add(seq.len(), Ordering::Relaxed);
                    }

                    if atty::is(atty::Stream::Stderr) {
                        eprint!(
                            "\rProcessed {} sequences ({} {})...",
                            processed_seq_ct.load(Ordering::Relaxed),
                            processed_ch_ct.load(Ordering::Relaxed),
                            if opts.input_is_protein { "aa" } else { "bp" }
                        );
                    }
                }
            });
        }
    });

    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }

    let final_seq_ct = processed_seq_ct.load(Ordering::Relaxed);
    let final_ch_ct = processed_ch_ct.load(Ordering::Relaxed);

    eprintln!(
        "Completed processing of {} sequences, {} {}",
        final_seq_ct,
        final_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );

    Ok(())
}

// Slightly slower but deterministic when multithreaded
fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, u64>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &mut Taxonomy,
) -> Result<()> {
    let mut processed_seq_ct = 0;
    let mut processed_ch_ct = 0;
    let mut sequence = Sequence::default();

    taxonomy.generate_external_to_internal_id_map();

    let stdin = io::stdin();
    let mut reader = BatchSequenceReader::new(stdin.lock());

    while reader.load_block(opts.block_size)? {
        while reader.next_sequence(&mut sequence)? {
            let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
            let mut taxid = 0;
            for seqid in &all_sequence_ids {
                if let Some(&ext_taxid) = id_to_taxon_map.get(seqid) {
                    if ext_taxid != 0 {
                        taxid = taxonomy
                            .lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                    }
                }
            }

            if taxid != 0 {
                let mut seq = sequence.seq.clone();
                if opts.input_is_protein && !seq.ends_with('*') {
                    seq.push('*');
                }

                process_sequence_deterministic(&seq, taxid, opts, kraken_index, taxonomy)?;

                processed_seq_ct += 1;
                processed_ch_ct += seq.len();

                if atty::is(atty::Stream::Stderr) {
                    eprint!(
                        "\rProcessed {} sequences ({} {})...",
                        processed_seq_ct,
                        processed_ch_ct,
                        if opts.input_is_protein { "aa" } else { "bp" }
                    );
                }
            }
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

// Process sequences fast in a single pass
fn process_sequence_fast_impl(
    seq: &str,
    taxid: u64,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
    scanner: &mut MinimizerScanner,
    min_clear_hash_value: u64,
) {
    scanner.load_sequence(seq, 0, seq.len());
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0 {
            let hc = murmur_hash3_64(&minimizer.to_le_bytes());
            if hc < min_clear_hash_value {
                continue;
            }
        }

        let mut existing_taxid = 0;
        let mut new_taxid = taxid;
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            new_taxid = tax.lowest_common_ancestor(new_taxid, existing_taxid);
        }
    }
}

// Process sequence deterministically to match C++ version
fn process_sequence_deterministic(
    seq: &str,
    taxid: u64,
    opts: &Options,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
) -> Result<()> {
    const SET_CT: usize = 256;
    let minimizer_sets = Arc::new(
        (0..SET_CT)
            .map(|_| Mutex::new(std::collections::BTreeSet::new()))
            .collect::<Vec<_>>(),
    );
    let hash = Arc::new(Mutex::new(hash));
    let tax = Arc::new(tax);

    // For each block in the sequence
    let mut j = 0;
    while j < seq.len() {
        let block_start = j;
        let mut block_finish = j + opts.block_size + opts.k - 1;
        if block_finish > seq.len() {
            block_finish = seq.len();
        }

        // Process each subblock in parallel
        let subblocks: Vec<_> = (block_start..block_finish)
            .step_by(opts.subblock_size)
            .collect();

        let minimizer_sets = Arc::clone(&minimizer_sets);
        let hash = Arc::clone(&hash);
        let tax = Arc::clone(&tax);

        subblocks.par_iter().for_each(|&i| {
            let mut scanner = MinimizerScanner::new(
                opts.k,
                opts.l,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
                true,
            );

            let mut subblock_finish = i + opts.subblock_size + opts.k - 1;
            if subblock_finish > block_finish {
                subblock_finish = block_finish;
            }

            scanner.load_sequence(seq, i, subblock_finish - i);

            while let Some(minimizer) = scanner.next_minimizer() {
                if scanner.is_ambiguous() {
                    continue;
                }

                let hash_code = murmur_hash3_64(&minimizer.to_le_bytes());
                if opts.min_clear_hash_value != 0 && hash_code < opts.min_clear_hash_value {
                    continue;
                }

                let zone = (hash_code % SET_CT as u64) as usize;
                if let Ok(mut set) = minimizer_sets[zone].lock() {
                    set.insert(minimizer);
                }
            }
        });

        // Combine all sets into sorted list
        let mut all_minimizers = Vec::new();
        let mut minimizer_lists = Vec::with_capacity(SET_CT);
        let mut prefix_sizes = vec![0; SET_CT + 1];

        // Convert sets to sorted vectors
        for (i, set_mutex) in minimizer_sets.iter().enumerate() {
            if let Ok(set) = set_mutex.lock() {
                let mut list: Vec<_> = set.iter().copied().collect();
                list.sort_unstable();
                prefix_sizes[i + 1] = list.len();
                minimizer_lists.push(list);
            }
            if let Ok(mut set) = set_mutex.lock() {
                set.clear(); // Clear set for next block
            }
        }

        // Calculate prefix sums
        for i in 2..=SET_CT {
            prefix_sizes[i] += prefix_sizes[i - 1];
        }

        // Initialize final list with correct size
        all_minimizers.resize(prefix_sizes[SET_CT], 0);

        // Copy sorted lists into final list
        for (i, list) in minimizer_lists.iter().enumerate() {
            let start = prefix_sizes[i];
            all_minimizers[start..start + list.len()].copy_from_slice(list);
        }

        // Process minimizers deterministically
        let mut mm_ct = all_minimizers.len();
        let mut index_list = Vec::with_capacity(mm_ct);
        let mut insertion_list = vec![true; mm_ct];

        // Loop to enforce deterministic order
        while !all_minimizers.is_empty() {
            // Gather insertion point information for all remaining minimizers
            for i in 0..mm_ct {
                if insertion_list[i] {
                    let mut idx = 0;
                    let hash_guard = hash.lock().unwrap();
                    insertion_list[i] = !hash_guard.find_index(all_minimizers[i], &mut idx);
                    index_list.push(idx);
                }
            }

            // Determine safe prefix to insert in parallel
            let mut novel_insertion_points = std::collections::BTreeSet::new();
            let mut safe_ct = 0;

            for i in 0..mm_ct {
                if insertion_list[i] {
                    if novel_insertion_points.contains(&index_list[i]) {
                        break;
                    }
                    novel_insertion_points.insert(index_list[i]);
                    safe_ct = i + 1;
                }
            }

            // Process safe prefix in parallel
            all_minimizers[..safe_ct].par_iter().for_each(|&minimizer| {
                if let Ok(mut hash_guard) = hash.lock() {
                    let mut existing_taxid = 0;
                    let mut new_taxid = taxid;
                    while !hash_guard.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
                        new_taxid = tax.lowest_common_ancestor(new_taxid, existing_taxid);
                    }
                }
            });

            // Remove processed prefix and re-iterate to process remainder
            all_minimizers = all_minimizers.split_off(safe_ct);
            index_list = index_list.split_off(safe_ct);
            insertion_list = insertion_list.split_off(safe_ct);
            mm_ct = all_minimizers.len();
        }

        j = block_finish - opts.k + 1;
    }

    Ok(())
}

fn main() -> Result<()> {
    let mut opts = parse_command_line()?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()?;

    let id_to_taxon_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;

    // Create NCBI taxonomy and convert to Kraken format
    eprint!("Creating taxonomy... ");
    let taxonomy = Taxonomy::from_ncbi_dmp(
        &format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory),
        &format!("{}/names.dmp", opts.ncbi_taxonomy_directory),
    )
    .context("Failed to create taxonomy from NCBI dump files")?;

    // Mark nodes used in the mapping
    for taxid in id_to_taxon_map.values() {
        if *taxid != 0 {
            NCBITaxonomyOps::mark_node(&taxonomy, *taxid);
        }
    }

    // Convert and save taxonomy
    NCBITaxonomyOps::convert_to_kraken_taxonomy(&taxonomy, &opts.taxonomy_filename)
        .context("Failed to convert taxonomy to Kraken format")?;

    eprintln!("done.");

    // Load the converted taxonomy
    eprint!("Loading taxonomy... ");
    let mut taxonomy = Taxonomy::from_kraken_taxonomy(&opts.taxonomy_filename)
        .context("Failed to load Kraken taxonomy")?;

    // Determine bits needed for the entire taxonomy
    let mut bits_needed_for_value = 1usize;
    while (1 << bits_needed_for_value) < taxonomy.node_count() {
        bits_needed_for_value += 1;
    }

    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        bail!(
            "More bits required for storing taxid than requested ({} > {})",
            bits_needed_for_value,
            opts.requested_bits_for_taxid
        );
    }
    let mut bits_for_taxid = bits_needed_for_value;
    if bits_for_taxid < opts.requested_bits_for_taxid {
        bits_for_taxid = opts.requested_bits_for_taxid;
    }

    // Handle MiniKraken build if needed
    let mut actual_capacity = opts.capacity;
    if opts.maximum_capacity > 0 {
        let frac = (opts.maximum_capacity as f64) / (opts.capacity as f64);
        if frac > 1.0 {
            bail!(
                "Maximum capacity {} larger than requested capacity {}",
                opts.maximum_capacity,
                opts.capacity
            );
        }
        opts.min_clear_hash_value = ((1.0 - frac) * (u64::MAX as f64)) as u64;
        actual_capacity = opts.maximum_capacity;
    }

    let mut kraken_index = CompactHashTable::new(
        actual_capacity,
        (32 - bits_for_taxid) as u8,
        bits_for_taxid as u8,
    );
    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    // Process sequences based on build mode
    if opts.deterministic_build {
        process_sequences(&opts, &id_to_taxon_map, &mut kraken_index, &mut taxonomy)?;
    } else {
        process_sequences_fast(&opts, &id_to_taxon_map, &mut kraken_index, &mut taxonomy)?;
    }

    // Write final data
    eprint!("Writing data to disk... ");
    kraken_index
        .write_table(&opts.hashtable_filename)
        .context("Failed to write hash table")?;

    let index_opts = IndexOptions {
        k: opts.k,
        l: opts.l,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        minimum_acceptable_hash_value: opts.min_clear_hash_value,
        revcom_version: CURRENT_REVCOM_VERSION,
        db_version: 0,
        db_type: 0,
    };

    let mut opts_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&opts.options_filename)
        .context("Failed to create options file")?;

    let bytes = unsafe {
        std::slice::from_raw_parts(
            &index_opts as *const IndexOptions as *const u8,
            std::mem::size_of::<IndexOptions>(),
        )
    };
    opts_file
        .write_all(bytes)
        .context("Failed to write options")?;
    eprintln!("complete.");

    Ok(())
}

fn expand_spaced_seed_mask(mask: &mut u64, bits_per_char: u8) {
    let mut new_mask = 0u64;
    let mut pos = 0;
    let mut old_pos = 0;

    while old_pos < 64 {
        if *mask & 1u64 << old_pos != 0 {
            new_mask |= 1u64 << pos;
            pos += bits_per_char as u64;
        }
        old_pos += 1;
    }
    *mask = new_mask;
}

pub fn set_minimizer_lca(hash: &mut CompactHashTable, minimizer: u64, taxid: u64, tax: &Taxonomy) {
    let mut old_value = 0;
    let mut new_value = taxid;
    while !hash.compare_and_set(minimizer, new_value, &mut old_value) {
        new_value = tax.lowest_common_ancestor(old_value, taxid);
    }
}

trait CompactHashOps {
    fn find_index(&self, key: u64, idx: &mut usize) -> bool;
    fn compare_and_set(&mut self, key: u64, new_value: u64, old_value: &mut u64) -> bool;
    fn write_table(&self, path: &str) -> Result<()>;
}

impl CompactHashOps for CompactHashTable {
    fn find_index(&self, key: u64, idx: &mut usize) -> bool {
        let probe_idx = (key & (self.node_count() as u64 - 1)) as usize;
        *idx = probe_idx;
        let value = self.get_at(probe_idx).load(Ordering::Relaxed);
        value != 0
    }

    fn compare_and_set(&mut self, key: u64, new_value: u64, old_value: &mut u64) -> bool {
        let probe_idx = (key & (self.node_count() as u64 - 1)) as usize;
        let current = self.get_at(probe_idx).load(Ordering::Relaxed);
        *old_value = current;
        if current == 0 {
            self.get_at(probe_idx).store(new_value, Ordering::Relaxed);
            return true;
        }
        false
    }

    fn write_table(&self, path: &str) -> Result<()> {
        let mut file = File::create(path)?;
        for i in 0..self.node_count() {
            let value = self.get_at(i).load(Ordering::Relaxed);
            file.write_all(&value.to_le_bytes())?;
        }
        Ok(())
    }
}
