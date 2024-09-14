/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::process;
use std::sync::{Arc, Mutex};

use clap::{App, Arg};
use rayon::prelude::*;

const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

#[derive(Debug)]
struct Options {
    id_to_taxon_map_filename: String,
    ncbi_taxonomy_directory: String,
    hashtable_filename: String,
    options_filename: String,
    taxonomy_filename: String,
    block_size: usize,
    subblock_size: usize,
    requested_bits_for_taxid: usize,
    num_threads: usize,
    input_is_protein: bool,
    k: usize,
    l: usize,
    capacity: usize,
    maximum_capacity: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
    deterministic_build: bool,
}

fn main() -> io::Result<()> {
    let mut opts = Options {
        spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
        toggle_mask: DEFAULT_TOGGLE_MASK,
        input_is_protein: false,
        num_threads: 1,
        block_size: DEFAULT_BLOCK_SIZE,
        subblock_size: DEFAULT_SUBBLOCK_SIZE,
        requested_bits_for_taxid: 0,
        min_clear_hash_value: 0,
        maximum_capacity: 0,
        deterministic_build: true,
        id_to_taxon_map_filename: String::new(),
        ncbi_taxonomy_directory: String::new(),
        hashtable_filename: String::new(),
        options_filename: String::new(),
        taxonomy_filename: String::new(),
        k: 0,
        l: 0,
        capacity: 0,
    };

    parse_command_line(&mut opts)?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    let id_to_taxon_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;

    generate_taxonomy(&opts, &id_to_taxon_map)?;

    eprintln!("Taxonomy parsed and converted.");

    let mut taxonomy = Taxonomy::new(&opts.taxonomy_filename, false)?;
    taxonomy.generate_external_to_internal_id_map();

    let mut bits_needed_for_value = 1;
    while (1 << bits_needed_for_value) < taxonomy.node_count() {
        bits_needed_for_value += 1;
    }
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        eprintln!("More bits required for storing taxid");
        process::exit(1);
    }

    let mut bits_for_taxid = bits_needed_for_value;
    if bits_for_taxid < opts.requested_bits_for_taxid {
        bits_for_taxid = opts.requested_bits_for_taxid;
    }

    let actual_capacity = if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            eprintln!("Maximum capacity larger than requested capacity");
            process::exit(1);
        }
        opts.min_clear_hash_value = ((1.0 - frac) * u64::MAX as f64) as u64;
        opts.maximum_capacity
    } else {
        opts.capacity
    };

    let mut kraken_index = CompactHashTable::new(
        actual_capacity,
        32 - bits_for_taxid as u32,
        bits_for_taxid as u32,
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
        k: opts.k as u32,
        l: opts.l as u32,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        minimum_acceptable_hash_value: opts.min_clear_hash_value,
        revcom_version: CURRENT_REVCOM_VERSION,
        db_version: 0,
        db_type: 0,
    };

    let mut opts_fs = File::create(&opts.options_filename)?;
    opts_fs.write_all(&bincode::serialize(&index_opts).unwrap())?;
    opts_fs.flush()?;

    eprintln!(" complete.");

    Ok(())
}

fn parse_command_line(opts: &mut Options) -> Result<(), io::Error> {
    let matches = App::new("build_db")
        .version("1.0")
        .about("Kraken 2 Database Builder")
        .arg(
            Arg::with_name("hashtable_filename")
                .short("H")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 hash table filename"),
        )
        .arg(
            Arg::with_name("id_to_taxon_map_filename")
                .short("m")
                .takes_value(true)
                .required(true)
                .help("Sequence ID to taxon map filename"),
        )
        .arg(
            Arg::with_name("taxonomy_filename")
                .short("t")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 taxonomy filename"),
        )
        .arg(
            Arg::with_name("ncbi_taxonomy_directory")
                .short("n")
                .takes_value(true)
                .required(true)
                .help("NCBI taxonomy directory name"),
        )
        .arg(
            Arg::with_name("options_filename")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 options filename"),
        )
        .arg(
            Arg::with_name("k")
                .short("k")
                .takes_value(true)
                .required(true)
                .help("Set length of k-mers"),
        )
        .arg(
            Arg::with_name("l")
                .short("l")
                .takes_value(true)
                .required(true)
                .help("Set length of minimizers"),
        )
        .arg(
            Arg::with_name("capacity")
                .short("c")
                .takes_value(true)
                .required(true)
                .help("Set capacity of hash table"),
        )
        .arg(
            Arg::with_name("maximum_capacity")
                .short("M")
                .takes_value(true)
                .help("Set maximum capacity of hash table (MiniKraken)"),
        )
        .arg(
            Arg::with_name("spaced_seed_mask")
                .short("S")
                .takes_value(true)
                .help("Spaced seed mask"),
        )
        .arg(
            Arg::with_name("toggle_mask")
                .short("T")
                .takes_value(true)
                .help("Minimizer toggle mask"),
        )
        .arg(
            Arg::with_name("input_is_protein")
                .short("X")
                .help("Input sequences are proteins"),
        )
        .arg(
            Arg::with_name("num_threads")
                .short("p")
                .takes_value(true)
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("deterministic_build")
                .short("F")
                .help("Use fast, nondeterministic building method"),
        )
        .arg(
            Arg::with_name("block_size")
                .short("B")
                .takes_value(true)
                .help("Read block size"),
        )
        .arg(
            Arg::with_name("subblock_size")
                .short("b")
                .takes_value(true)
                .help("Read subblock size"),
        )
        .arg(
            Arg::with_name("requested_bits_for_taxid")
                .short("r")
                .takes_value(true)
                .help("Bit storage requested for taxid"),
        )
        .get_matches();

    opts.hashtable_filename = matches.value_of("hashtable_filename").unwrap().to_string();
    opts.id_to_taxon_map_filename = matches
        .value_of("id_to_taxon_map_filename")
        .unwrap()
        .to_string();
    opts.taxonomy_filename = matches.value_of("taxonomy_filename").unwrap().to_string();
    opts.ncbi_taxonomy_directory = matches
        .value_of("ncbi_taxonomy_directory")
        .unwrap()
        .to_string();
    opts.options_filename = matches.value_of("options_filename").unwrap().to_string();
    opts.k = matches.value_of("k").unwrap().parse().unwrap();
    opts.l = matches.value_of("l").unwrap().parse().unwrap();
    opts.capacity = matches.value_of("capacity").unwrap().parse().unwrap();

    if let Some(value) = matches.value_of("maximum_capacity") {
        opts.maximum_capacity = value.parse().unwrap();
    }

    if let Some(value) = matches.value_of("spaced_seed_mask") {
        opts.spaced_seed_mask = u64::from_str_radix(value, 2).unwrap();
    }

    if let Some(value) = matches.value_of("toggle_mask") {
        opts.toggle_mask = u64::from_str_radix(value, 2).unwrap();
    }

    if matches.is_present("input_is_protein") {
        opts.input_is_protein = true;
    }

    if let Some(value) = matches.value_of("num_threads") {
        opts.num_threads = value.parse().unwrap();
    }

    if matches.is_present("deterministic_build") {
        opts.deterministic_build = false;
    }

    if let Some(value) = matches.value_of("block_size") {
        opts.block_size = value.parse().unwrap();
    }

    if let Some(value) = matches.value_of("subblock_size") {
        opts.subblock_size = value.parse().unwrap();
    }

    if let Some(value) = matches.value_of("requested_bits_for_taxid") {
        opts.requested_bits_for_taxid = value.parse().unwrap();
    }

    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        process::exit(1);
    }

    if opts.block_size < opts.subblock_size {
        eprintln!("block size cannot be less than subblock size");
        process::exit(1);
    }

    if opts.maximum_capacity > opts.capacity {
        eprintln!("maximum capacity option shouldn't specify larger capacity than normal");
        process::exit(1);
    }

    Ok(())
}

fn read_id_to_taxon_map(filename: &str) -> io::Result<HashMap<String, TaxId>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            let taxid: TaxId = taxid_str.parse().unwrap();
            if taxid != 0 {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }

    Ok(id_map)
}

fn generate_taxonomy(opts: &Options, id_map: &HashMap<String, TaxId>) -> io::Result<()> {
    let nodes_filename = format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory);
    let names_filename = format!("{}/names.dmp", opts.ncbi_taxonomy_directory);

    let mut ncbi_taxonomy = NCBITaxonomy::new(&nodes_filename, &names_filename)?;

    for &taxid in id_map.values() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid);
        }
    }

    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename)
}

fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, TaxId>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &Taxonomy,
) -> io::Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    let stdin = io::stdin();
    let mut reader = stdin.lock();

    let mut sequence_reader = BatchSequenceReader::new();
    let mut sequence = Sequence::default();

    while sequence_reader.load_block(&mut reader, opts.block_size)? {
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
                let seq = if opts.input_is_protein && !sequence.seq.ends_with('*') {
                    format!("{}*", sequence.seq)
                } else {
                    sequence.seq.clone()
                };

                process_sequence(opts, &seq, taxid, kraken_index, taxonomy, &sequence);

                processed_seq_ct += 1;
                processed_ch_ct += sequence.seq.len();
            }
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

fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    let chars: Vec<char> = header.chars().collect();
    let mut i = 1; // Start after '>'
    while i < chars.len() {
        if chars[i] == '\x01' {
            // 0x01 starts new ID at next char
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = true;
        } else if in_id && chars[i].is_whitespace() {
            // Spaces end ID
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

fn process_sequence(
    opts: &Options,
    seq: &str,
    taxid: TaxId,
    kraken_index: &mut CompactHashTable,
    taxonomy: &Taxonomy,
    sequence: &Sequence,
) {
    let set_ct = 256;
    let mut locks = Vec::with_capacity(set_ct);
    for _ in 0..set_ct {
        locks.push(Mutex::new(()));
    }

    let seq_len = seq.len();
    let mut blocks = Vec::new();

    let k = opts.k;
    let block_size = opts.block_size;
    let subblock_size = opts.subblock_size;

    // Divide sequence into blocks
    let mut j = 0;
    while j < seq_len {
        let block_start = j;
        let mut block_finish = j + block_size + k - 1;
        if block_finish > seq_len {
            block_finish = seq_len;
        }
        blocks.push((block_start, block_finish));
        j += block_size;
    }

    for (block_start, block_finish) in blocks {
        let mut minimizer_sets: Vec<HashSet<u64>> = vec![HashSet::new(); set_ct];

        // Process subblocks in parallel
        (block_start..block_finish)
            .into_par_iter()
            .step_by(subblock_size)
            .for_each(|i| {
                let subblock_finish = if i + subblock_size + k - 1 > block_finish {
                    block_finish
                } else {
                    i + subblock_size + k - 1
                };

                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    opts.spaced_seed_mask,
                    !opts.input_is_protein,
                    opts.toggle_mask,
                    CURRENT_REVCOM_VERSION,
                );

                scanner.load_sequence(&seq[i..subblock_finish]);

                while let Some(minimizer) = scanner.next_minimizer() {
                    if scanner.is_ambiguous() {
                        continue;
                    }
                    let hc = murmur_hash3(minimizer);
                    if opts.min_clear_hash_value != 0 && hc < opts.min_clear_hash_value {
                        continue;
                    }
                    let zone = (hc % set_ct as u64) as usize;
                    let _lock = locks[zone].lock().unwrap();
                    minimizer_sets[zone].insert(minimizer);
                }
            });

        // Combine sets into sorted list
        let mut minimizer_list = Vec::new();
        for minimizer_set in minimizer_sets {
            let mut minimizers: Vec<u64> = minimizer_set.into_iter().collect();
            minimizers.sort_unstable();
            minimizer_list.extend(minimizers);
        }

        // Remove duplicates
        minimizer_list.sort_unstable();
        minimizer_list.dedup();

        // Process minimizers
        for &minimizer in &minimizer_list {
            set_minimizer_lca(kraken_index, minimizer, taxid, taxonomy);
        }
    }
}

fn set_minimizer_lca(
    hash: &mut CompactHashTable,
    minimizer: u64,
    taxid: TaxId,
    taxonomy: &Taxonomy,
) {
    let mut existing_taxid = 0;
    let mut new_taxid = taxid;

    while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
        new_taxid = taxonomy.lowest_common_ancestor(new_taxid, existing_taxid);
    }
}

fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, TaxId>,
    kraken_index: &mut CompactHashTable,
    taxonomy: &Taxonomy,
) -> io::Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    let stdin = io::stdin();
    let reader = stdin.lock();

    let id_to_taxon_map = Arc::new(id_to_taxon_map.clone());
    let kraken_index = Arc::new(Mutex::new(kraken_index));
    let taxonomy = Arc::new(taxonomy.clone());

    let sequences: Vec<Sequence> = read_sequences(reader)?;

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

            scanner.load_sequence(&seq);

            while let Some(minimizer) = scanner.next_minimizer() {
                if scanner.is_ambiguous() {
                    continue;
                }
                if opts.min_clear_hash_value != 0
                    && murmur_hash3(minimizer) < opts.min_clear_hash_value
                {
                    continue;
                }
                let mut kraken_index = kraken_index.lock().unwrap();
                set_minimizer_lca(&mut kraken_index, minimizer, taxid, &taxonomy);
            }

            processed_seq_ct += 1;
            processed_ch_ct += seq.len();
        }
    });

    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );

    Ok(())
}

fn read_sequences<R: BufRead>(reader: R) -> io::Result<Vec<Sequence>> {
    let mut sequences = Vec::new();
    let mut sequence = Sequence::default();
    let mut lines = reader.lines();

    while let Some(line) = lines.next() {
        let line = line?;
        if line.starts_with('>') {
            if !sequence.seq.is_empty() {
                sequences.push(sequence.clone());
                sequence = Sequence::default();
            }
            sequence.header = line;
        } else {
            sequence.seq.push_str(&line);
        }
    }

    if !sequence.seq.is_empty() {
        sequences.push(sequence);
    }

    Ok(sequences)
}
