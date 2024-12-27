//
// Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
//
// This file is part of the Kraken 2 taxonomic sequence classification system.
//

//! This Rust module is the translated version of `build_db.cc` from Kraken 2.
//! It handles command-line parsing, reading sequences, building a compact hash table,
//! and finalizing the Kraken 2 database building process.

use std::collections::{BTreeSet, HashMap};
use std::env;
use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::exit;
use std::sync::{Arc, Mutex};

use anyhow::{anyhow, bail, Context, Result};
use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;

pub const BITS_PER_CHAR_DNA: u8 = 2;
pub const BITS_PER_CHAR_PRO: u8 = 5;

// Dummy modules that represent functionality analogous to the C++ includes.
// In a real-world scenario, you would implement each module's logic, or
// replace them with your own crate/module references.
use crate::compact_hash::{set_minimizer_lca, CompactHashTable};
use crate::kraken2_data::{
    IndexOptions, CURRENT_REVCOM_VERSION, DEFAULT_SPACED_SEED_MASK, DEFAULT_TOGGLE_MASK,
};

use crate::mmscanner::{murmur_hash3_64, MinimizerScanner};
use crate::seqreader::{BatchSequenceReader, Sequence};
use crate::taxonomy::Taxonomy;
use crate::utilities::expand_spaced_seed_mask as ExpandSpacedSeedMask;

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
    pub k: isize,
    pub l: isize,
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
        .map_or(0, |s| s.parse::<isize>().unwrap_or(0));
    let l = matches
        .get_one::<String>("l")
        .map_or(0, |s| s.parse::<isize>().unwrap_or(0));
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
    let spaced_seed_mask = if let Some(s) = spaced_seed_mask_str {
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

    if k < l as isize {
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
        ExpandSpacedSeedMask(spaced_seed_mask, bits_per_char);
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
    // This function deals with NCBI's usage of 0x01 to denote
    // multiple IDs in one FASTA header
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;
    // Start from character after '>' in typical FASTA header
    for (i, c) in header.chars().enumerate().skip(1) {
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

fn read_id_to_taxon_map(filename: &str) -> Result<HashMap<String, u32>> {
    let file =
        File::open(filename).with_context(|| format!("unable to read from '{}'", filename))?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            let taxid = taxid_str.parse::<u32>().unwrap_or(0);
            if taxid != 0 {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
    Ok(id_map)
}

fn generate_taxonomy(opts: &mut Options, id_map: &HashMap<String, u32>) -> Result<()> {
    // Build NCBI taxonomy from nodes.dmp and names.dmp, then convert to Kraken
    let ncbi_taxonomy = Taxonomy::from_ncbi_dmp(
        &format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory),
        &format!("{}/names.dmp", opts.ncbi_taxonomy_directory),
    )?;

    for taxid in id_map.values() {
        if *taxid != 0 {
            ncbi_taxonomy.mark_node(*taxid);
        }
    }
    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename)?;
    Ok(())
}

// A quick but nondeterministic build
fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, u32>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &Taxonomy,
) -> Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    // In Rust, we'll use rayon or standard threads instead of OpenMP
    // We'll read everything from STDIN in blocks
    // This is a simplified approach: we read line by line and process in parallel.

    let stdin = io::stdin();
    let mut reader = crate::seqreader::BatchSequenceReader::new(stdin.lock());

    // We'll emulate block reading with a while loop:
    loop {
        // Attempt to load next sequence in a single-thread manner
        let maybe_seq = reader.next_sequence()?;
        if maybe_seq.is_none() {
            break; // no more data
        }
        let sequence = maybe_seq.unwrap();

        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let mut taxid = 0u32;
        for seqid in all_sequence_ids {
            if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                if ext_taxid != 0 {
                    taxid =
                        taxonomy.lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                }
            }
        }
        if taxid != 0 {
            let mut seq_str = sequence.seq;
            if opts.input_is_protein && !seq_str.ends_with('*') {
                seq_str.push('*');
            }
            let mut scanner = MinimizerScanner::new(
                opts.k as usize,
                opts.l as usize,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
            );
            process_sequence_fast_impl(
                &seq_str,
                taxid,
                kraken_index,
                taxonomy,
                &mut scanner,
                opts.min_clear_hash_value,
            );
            processed_seq_ct += 1;
            processed_ch_ct += seq_str.len();
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

// Slightly slower but deterministic when multithreaded
fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, u32>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &Taxonomy,
) -> Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    let stdin = io::stdin();
    let mut reader = super::seqreader::BatchSequenceReader::new(stdin.lock());

    while let Some(sequence) = reader.next_sequence()? {
        let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
        let mut taxid = 0u32;
        for seqid in all_sequence_ids {
            if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                if ext_taxid != 0 {
                    taxid =
                        taxonomy.lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                }
            }
        }
        if taxid != 0 {
            let mut seq_str = sequence.seq;
            if opts.input_is_protein && !seq_str.ends_with('*') {
                seq_str.push('*');
            }
            process_sequence(opts, &seq_str, taxid, kraken_index, taxonomy)?;
            processed_seq_ct += 1;
            processed_ch_ct += seq_str.len();
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

// This function exists to mimic the C++ version which sets LCA in a single pass:
fn process_sequence_fast_impl(
    seq: &str,
    taxid: u32,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
    scanner: &mut MinimizerScanner,
    min_clear_hash_value: u64,
) {
    scanner.load_sequence(seq);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0
            && murmur_hash3_64(&minimizer.to_le_bytes()) < min_clear_hash_value
        {
            continue;
        }
        let mut existing_taxid = 0;
        let mut new_taxid = taxid;
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            new_taxid = tax.lowest_common_ancestor(new_taxid, existing_taxid);
        }
    }
}

// This function parallels the C++ `ProcessSequence()`, ensuring a deterministic insertion order.
fn process_sequence(
    opts: &Options,
    seq: &str,
    taxid: u32,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
) -> Result<()> {
    let set_ct = 256;
    let mut locks = Vec::with_capacity(set_ct);
    for _ in 0..set_ct {
        locks.push(Mutex::new(()));
    }
    // for each block in the sequence
    let bytes = seq.as_bytes();

    let mut j = 0;
    while j < bytes.len() {
        let block_start = j;
        let mut block_finish = j + opts.block_size + (opts.k as usize) - 1;
        if block_finish > bytes.len() {
            block_finish = bytes.len();
        }
        let mut minimizer_sets: Vec<Mutex<BTreeSet<u64>>> = Vec::with_capacity(set_ct);
        for _ in 0..set_ct {
            minimizer_sets.push(Mutex::new(BTreeSet::new()));
        }

        // Dividing into subblocks in parallel in C++ was done via OpenMP;
        // here we can do something with rayon, but to keep it simpler, we do a normal loop
        let mut i = block_start;
        while i < block_finish {
            let subblock_finish =
                std::cmp::min(i + opts.subblock_size + (opts.k as usize) - 1, block_finish);
            let mut scanner = MinimizerScanner::new(
                opts.k as usize,
                opts.l as usize,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
            );
            scanner.load_sequence(&bytes[i..subblock_finish]);
            while let Some(minimizer) = scanner.next_minimizer() {
                if scanner.is_ambiguous() {
                    continue;
                }
                let hc = murmur_hash3_64(&minimizer.to_le_bytes());
                if opts.min_clear_hash_value != 0 && hc < opts.min_clear_hash_value {
                    continue;
                }
                let zone = (hc % set_ct as u64) as usize;
                {
                    let _lk = locks[zone].lock().unwrap();
                    let mut set_ref = minimizer_sets[zone].lock().unwrap();
                    set_ref.insert(minimizer);
                }
            }
            i += opts.subblock_size;
        }

        // combine sets into one sorted list
        let mut combined = Vec::new();
        for zone in 0..set_ct {
            let set_guard = minimizer_sets[zone].lock().unwrap();
            combined.extend(set_guard.iter().copied());
        }

        // In place sort
        combined.sort_unstable();

        // We'll store parallel arrays to track insertion
        let mm_ct = combined.len();
        let mut insertion_list = vec![true; mm_ct];
        let mut index_list = vec![0usize; mm_ct];

        // Deterministic insertion loop
        let mut idx = 0;
        while idx < combined.len() {
            // Gather insertion point information
            for i in idx..combined.len() {
                if insertion_list[i] {
                    let mut found_idx = 0usize;
                    let found = hash.find_index(combined[i], &mut found_idx);
                    if found {
                        // Already present in the table
                        insertion_list[i] = false;
                    }
                    index_list[i] = found_idx;
                }
            }
            // Determine safe prefix
            let mut novel_insertion_points = BTreeSet::new();
            let mut safe_ct = 0;
            for i in idx..combined.len() {
                if insertion_list[i] {
                    if novel_insertion_points.contains(&index_list[i]) {
                        break;
                    }
                    novel_insertion_points.insert(index_list[i]);
                    safe_ct += 1;
                } else {
                    safe_ct += 1;
                }
            }

            // Actually insert them
            let end_pos = idx + safe_ct;
            for i in idx..end_pos {
                if insertion_list[i] {
                    set_minimizer_lca(hash, combined[i], taxid, tax);
                }
            }
            idx += safe_ct;
        }

        j = block_finish - (opts.k as usize) + 1;
    }

    Ok(())
}

fn main() -> Result<()> {
    let mut opts = parse_command_line()?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()?;

    let id_to_taxon_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;
    generate_taxonomy(&mut opts, &id_to_taxon_map)?;

    eprintln!("Taxonomy parsed and converted.");

    let taxonomy = Taxonomy::from_kraken_taxonomy(&opts.taxonomy_filename)?;
    taxonomy.generate_external_to_internal_id_map();

    // Determine bits needed for the entire taxonomy
    let mut bits_needed_for_value = 1usize;
    while (1 << bits_needed_for_value) < taxonomy.node_count() {
        bits_needed_for_value += 1;
    }

    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        bail!("more bits required for storing taxid");
    }
    let mut bits_for_taxid = bits_needed_for_value;
    if bits_for_taxid < opts.requested_bits_for_taxid {
        bits_for_taxid = opts.requested_bits_for_taxid;
    }

    // Possibly build a MiniKraken
    let mut actual_capacity = opts.capacity;
    if opts.maximum_capacity > 0 {
        let frac = (opts.maximum_capacity as f64) / (opts.capacity as f64);
        if frac > 1.0 {
            bail!("maximum capacity larger than requested capacity");
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

    if opts.deterministic_build {
        process_sequences(&opts, &id_to_taxon_map, &mut kraken_index, &taxonomy)?;
    } else {
        process_sequences_fast(&opts, &id_to_taxon_map, &mut kraken_index, &taxonomy)?;
    }

    eprint!("Writing data to disk... ");
    kraken_index.write_table(&opts.hashtable_filename)?;

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

    // We'll write the IndexOptions struct in little-endian byte form
    // In a real scenario, you'd carefully replicate the same layout as the C++ code
    let mut opts_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&opts.options_filename)
        .context("Unable to create options file")?;
    let bytes = unsafe {
        std::slice::from_raw_parts(
            &index_opts as *const IndexOptions as *const u8,
            std::mem::size_of::<IndexOptions>(),
        )
    };
    opts_file.write_all(bytes)?;
    eprintln!("complete.");

    Ok(())
}
