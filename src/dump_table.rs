/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::HashMap;
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, Read, Write};
use std::process;

use clap::{App, Arg};

mod compact_hash;
mod kraken2_data;
mod kraken2_headers;
mod mmscanner;
mod reports;
mod taxonomy;
mod utilities;

use compact_hash::CompactHashTable;
use kraken2_data::*;
use kraken2_headers::*;
use reports::{report_kraken_style, report_mpa_style};
use taxonomy::Taxonomy;
use utilities::*;

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

fn main() -> io::Result<()> {
    let mut opts = Options {
        hashtable_filename: String::new(),
        taxonomy_filename: String::new(),
        options_filename: String::new(),
        output_filename: "/dev/stdout".to_string(),
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

    let kraken_index = CompactHashTable::new_from_file(&opts.hashtable_filename)?;
    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, false)?;

    let idx_opts = read_index_options(&opts.options_filename)?;

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
                * if idx_opts.dna_db {
                    BITS_PER_CHAR_DNA
                } else {
                    BITS_PER_CHAR_PRO
                }
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

        let taxid_counters: HashMap<TaxId, ReadCounter> = taxid_counts
            .iter()
            .map(|(&taxid, &count)| (taxid, ReadCounter { count, paired: 0 }))
            .collect();

        if opts.use_mpa_style {
            report_mpa_style(
                &opts.output_filename,
                opts.report_zeros,
                &taxonomy,
                &taxid_counters,
            )?;
        } else {
            report_kraken_style(
                &opts.output_filename,
                opts.report_zeros,
                false,
                &taxonomy,
                &taxid_counters,
                total_seqs,
                0,
            )?;
        }
    }

    Ok(())
}

fn parse_command_line(opts: &mut Options) -> Result<(), io::Error> {
    let matches = App::new("dump_table")
        .version("1.0")
        .about("Kraken 2 Database Dumper")
        .arg(
            Arg::with_name("hashtable_filename")
                .short("H")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 hash table filename"),
        )
        .arg(
            Arg::with_name("taxonomy_filename")
                .short("t")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 taxonomy filename"),
        )
        .arg(
            Arg::with_name("options_filename")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("Kraken 2 database options filename"),
        )
        .arg(
            Arg::with_name("output_filename")
                .short("O")
                .takes_value(true)
                .help("Output filename (default: /dev/stdout)"),
        )
        .arg(
            Arg::with_name("use_mpa_style")
                .short("m")
                .help("Use MPA style output instead of Kraken 2 output"),
        )
        .arg(
            Arg::with_name("skip_counts")
                .short("s")
                .help("Skip reporting minimizer counts, just show DB stats"),
        )
        .arg(
            Arg::with_name("report_zeros")
                .short("z")
                .help("Report taxa with zero counts"),
        )
        .arg(
            Arg::with_name("num_threads")
                .short("p")
                .takes_value(true)
                .help("Number of threads"),
        )
        .get_matches();

    opts.hashtable_filename = matches.value_of("hashtable_filename").unwrap().to_string();
    opts.taxonomy_filename = matches.value_of("taxonomy_filename").unwrap().to_string();
    opts.options_filename = matches.value_of("options_filename").unwrap().to_string();

    if let Some(value) = matches.value_of("output_filename") {
        opts.output_filename = value.to_string();
    }

    if matches.is_present("use_mpa_style") {
        opts.use_mpa_style = true;
    }

    if matches.is_present("skip_counts") {
        opts.skip_counts = true;
    }

    if matches.is_present("report_zeros") {
        opts.report_zeros = true;
    }

    if let Some(value) = matches.value_of("num_threads") {
        opts.num_threads = value.parse().unwrap();
    }

    Ok(())
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
