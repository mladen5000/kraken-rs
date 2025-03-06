/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};

use anyhow::Result;
use bincode;
use clap::{Arg, Command};
use rayon;

use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{IndexOptions, TaxId, BITS_PER_CHAR_DNA, BITS_PER_CHAR_PRO};
use crate::readcounts::{HyperLogLogPlusMinus, ReadCounts};
use crate::reports::{report_kraken_style, report_mpa_style};
use crate::taxonomy::Taxonomy;

const CURRENT_REVCOM_VERSION: u32 = 1;
const DEFAULT_BITS_FOR_TAXID: u8 = 32;

type TaxidCounters = HashMap<TaxId, ReadCounts<HyperLogLogPlusMinus>>;

#[derive(Debug)]
struct Options {
    hashtable_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    output_filename: String,
    use_mpa_style: bool,
    report_zeros: bool,
    skip_counts: bool,
    num_threads: usize,
}

fn main() -> Result<()> {
    let mut opts = Options {
        hashtable_filename: String::new(),
        taxonomy_filename: String::new(),
        options_filename: String::new(),
        output_filename: "/dev/fd/1".to_string(), // Changed to match C++ default
        use_mpa_style: false,
        report_zeros: false,
        skip_counts: false,
        num_threads: 1,
    };

    parse_command_line(&mut opts)?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    // Read the index options first to get value and taxid bits
    let idx_opts = read_index_options(&opts.options_filename)?;

    // Create and load the hash table
    let mut file = File::open(&opts.hashtable_filename)?;
    let metadata = file.metadata()?;
    let capacity = metadata.len() as usize / std::mem::size_of::<u64>();

    let bits_for_taxid = DEFAULT_BITS_FOR_TAXID;
    let value_bits = (32 - bits_for_taxid) as u8;
    let kraken_index = CompactHashTable::new(capacity, value_bits, bits_for_taxid);

    let mut buffer = vec![0u64; capacity];
    file.read_exact(unsafe {
        std::slice::from_raw_parts_mut(buffer.as_mut_ptr() as *mut u8, capacity * 8)
    })?;
    for (i, &value) in buffer.iter().enumerate() {
        kraken_index
            .get_at(i)
            .store(value, std::sync::atomic::Ordering::Relaxed);
    }

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, false)?;

    println!(
        "# Database options: {} db, k = {}, l = {}",
        if idx_opts.dna_db {
            "nucleotide"
        } else {
            "protein"
        },
        idx_opts.k,
        idx_opts.l
    );
    println!(
        "# Spaced mask = {}",
        mask_to_str(
            idx_opts.spaced_seed_mask,
            idx_opts.l as usize
                * (if idx_opts.dna_db {
                    BITS_PER_CHAR_DNA as usize
                } else {
                    BITS_PER_CHAR_PRO as usize
                })
        )
    );
    println!("# Toggle mask = {}", mask_to_str(idx_opts.toggle_mask, 64));
    println!("# Total taxonomy nodes: {}", taxonomy.node_count());
    println!("# Table size: {}", kraken_index.size());
    println!("# Table capacity: {}", kraken_index.capacity());
    println!(
        "# Min clear hash value = {}",
        idx_opts.minimum_acceptable_hash_value
    );
    if idx_opts.revcom_version != CURRENT_REVCOM_VERSION {
        println!("# Built with outdated revcom version");
    }

    if !opts.skip_counts {
        let taxid_counts = kraken_index.get_value_counts();
        let total_seqs: u64 = taxid_counts.values().sum();

        let taxid_counters: TaxidCounters = taxid_counts
            .iter()
            .map(|(&taxid, &count)| (taxid, ReadCounts::with_counts(count, count)))
            .collect();

        if opts.use_mpa_style {
            report_mpa_style(
                &opts.output_filename,
                opts.report_zeros,
                &taxonomy,
                &taxid_counters,
            )
            .map_err(|e| anyhow::anyhow!("Failed to generate MPA report: {}", e))?;
        } else {
            report_kraken_style(
                &opts.output_filename,
                opts.report_zeros,
                false,
                &taxonomy,
                &taxid_counters,
                total_seqs,
                0,
            )
            .map_err(|e| anyhow::anyhow!("Failed to generate Kraken report: {}", e))?;
        }
    }

    Ok(())
}

fn parse_command_line(opts: &mut Options) -> Result<(), io::Error> {
    let matches = Command::new("dump_table")
        .version("1.0")
        .about("Kraken 2 Database Dumper")
        .arg(
            Arg::new("hashtable_filename")
                .short('H')
                .required(true)
                .value_parser(clap::value_parser!(String))
                .help("Kraken 2 hash table filename"),
        )
        .arg(
            Arg::new("taxonomy_filename")
                .short('t')
                .required(true)
                .value_parser(clap::value_parser!(String))
                .help("Kraken 2 taxonomy filename"),
        )
        .arg(
            Arg::new("options_filename")
                .short('o')
                .required(true)
                .value_parser(clap::value_parser!(String))
                .help("Kraken 2 database options filename"),
        )
        .arg(
            Arg::new("output_filename")
                .short('O')
                .value_parser(clap::value_parser!(String))
                .help("Output filename (def: /dev/fd/1)"),
        )
        .arg(
            Arg::new("use_mpa_style")
                .short('m')
                .action(clap::ArgAction::SetTrue)
                .help("Use MPA style output instead of Kraken 2 output"),
        )
        .arg(
            Arg::new("skip_counts")
                .short('s')
                .action(clap::ArgAction::SetTrue)
                .help("Skip reporting minimizer counts, just show DB stats"),
        )
        .arg(
            Arg::new("report_zeros")
                .short('z')
                .action(clap::ArgAction::SetTrue)
                .help("Report taxa with zero counts"),
        )
        .arg(
            Arg::new("num_threads")
                .short('p')
                .value_parser(clap::value_parser!(usize))
                .help("Number of threads"),
        )
        .get_matches();

    // Check required parameters
    opts.hashtable_filename = matches
        .get_one::<String>("hashtable_filename")
        .expect("hashtable_filename is required")
        .clone();

    opts.taxonomy_filename = matches
        .get_one::<String>("taxonomy_filename")
        .expect("taxonomy_filename is required")
        .clone();

    opts.options_filename = matches
        .get_one::<String>("options_filename")
        .expect("options_filename is required")
        .clone();

    if let Some(value) = matches.get_one::<String>("output_filename") {
        opts.output_filename = value.clone();
    }

    opts.use_mpa_style = matches.get_flag("use_mpa_style");
    opts.skip_counts = matches.get_flag("skip_counts");
    opts.report_zeros = matches.get_flag("report_zeros");

    if let Some(value) = matches.get_one::<usize>("num_threads") {
        opts.num_threads = *value;
    }

    Ok(())
}

fn print_usage() {
    eprintln!(
        "Usage: dump_table <options>\n\n\
        Options (*mandatory):\n\
        * -H FILENAME   Kraken 2 hash table filename\n\
        * -t FILENAME   Kraken 2 taxonomy filename\n\
        * -o FILENAME   Kraken 2 database options filename\n\
          -O FILENAME   Output filename (def: /dev/fd/1)\n\
          -m            Use MPA style output instead of Kraken 2 output\n\
          -s            Skip reporting minimizer counts, just show DB stats\n\
          -z            Report taxa with zero counts"
    );
}

fn read_index_options(filename: &str) -> io::Result<IndexOptions> {
    let mut file = File::open(filename)?;
    let mut buffer = vec![0; std::mem::size_of::<IndexOptions>()];
    file.read_exact(&mut buffer)?;
    let idx_opts: IndexOptions = bincode::deserialize(&buffer).unwrap();
    Ok(idx_opts)
}

fn mask_to_str(mask: u64, digits: usize) -> String {
    let mut str = String::new();
    for i in (0..digits).rev() {
        str.push(if (mask >> i) & 1 == 1 { '1' } else { '0' });
    }
    str
}
