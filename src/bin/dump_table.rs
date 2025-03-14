// Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
// Rust conversion Copyright 2025
//
// This file is part of the Kraken 2 taxonomic sequence classification system.

use std::fs::File;
use std::io::{self, Read, Write};
use std::process::exit;

use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;

use kraken_rs::compact_hash::CompactHashTable;
use kraken_rs::kraken2_data::{
    IndexOptions, TaxonCounters, TaxonCountersMap, BITS_PER_CHAR_DNA, BITS_PER_CHAR_PRO,
};
use kraken_rs::reports::{report_kraken_style, report_mpa_style};
use kraken_rs::taxonomy::Taxonomy;
use kraken_rs::utilities::CURRENT_REVCOM_VERSION;

// Define command-line options struct to match the C++ implementation
struct Options {
    hashtable_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    output_filename: String,
    use_mpa_style: bool,
    report_zeros: bool,
    skip_counts: bool,
    memory_mapping: bool,
    num_threads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            hashtable_filename: String::new(),
            taxonomy_filename: String::new(),
            options_filename: String::new(),
            output_filename: "/dev/fd/1".to_string(),
            use_mpa_style: false,
            report_zeros: false,
            skip_counts: false,
            memory_mapping: false,
            num_threads: 1,
        }
    }
}

// Convert a mask to a binary string
fn mask2str(mask: u64, digits: u32) -> String {
    (0..digits)
        .rev()
        .map(|i| if (mask >> i) & 1 == 1 { '1' } else { '0' })
        .collect()
}

// Parse command-line arguments
fn parse_command_line(args: &[String], opts: &mut Options) -> Result<()> {
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => {
                usage(0);
            }
            "-H" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("missing argument for -H");
                    usage(1);
                }
                opts.hashtable_filename = args[i].clone();
            }
            "-t" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("missing argument for -t");
                    usage(1);
                }
                opts.taxonomy_filename = args[i].clone();
            }
            "-o" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("missing argument for -o");
                    usage(1);
                }
                opts.options_filename = args[i].clone();
            }
            "-z" => {
                opts.report_zeros = true;
            }
            "-O" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("missing argument for -O");
                    usage(1);
                }
                opts.output_filename = args[i].clone();
            }
            "-m" => {
                opts.use_mpa_style = true;
            }
            "-s" => {
                opts.skip_counts = true;
            }
            "-p" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("missing argument for -p");
                    usage(1);
                }
                opts.num_threads = args[i].parse().unwrap_or(1);
            }
            "-M" => {
                opts.memory_mapping = true;
            }
            _ => {
                eprintln!("unknown option: {}", args[i]);
                usage(1);
            }
        }
        i += 1;
    }

    if opts.hashtable_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        eprintln!("missing mandatory filename parameter");
        usage(1);
    }

    Ok(())
}

// Print usage information
fn usage(exit_code: i32) -> ! {
    eprintln!(
        "Usage: dump_table <options>\n\n\
        Options (*mandatory):\n\
        * -H FILENAME   Kraken 2 hash table filename\n\
        * -t FILENAME   Kraken 2 taxonomy filename\n\
        * -o FILENAME   Kraken 2 database options filename\n\
        -O FILENAME   Output filename (def: /dev/fd/1)\n\
        -m            Use MPA style output instead of Kraken 2 output\n\
        -M            Use memory mapping to access hash & taxonomy\n\
        -s            Skip reporting minimizer counts, just show DB stats\n\
        -p INT        Number of threads\n\
        -z            Report taxa with zero counts"
    );
    exit(exit_code);
}

pub fn main() -> Result<()> {
    // Initialize options with defaults
    let mut opts = Options::default();

    // Parse command-line arguments
    let args: Vec<String> = std::env::args().collect();
    parse_command_line(&args, &mut opts)?;

    // Set number of threads
    ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    // Initialize the Kraken index and taxonomy
    let kraken_index = CompactHashTable::from_file(&opts.hashtable_filename, opts.memory_mapping)
        .context("Failed to load hash table")?;
    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.memory_mapping)
        .context("Failed to load taxonomy")?;

    // Read the index options
    let mut idx_opts = IndexOptions::default();
    let mut file = File::open(&opts.options_filename).context("Failed to open options file")?;

    // Read the binary option data
    unsafe {
        let idx_opts_ptr = &mut idx_opts as *mut _ as *mut u8;
        let idx_opts_slice =
            std::slice::from_raw_parts_mut(idx_opts_ptr, std::mem::size_of::<IndexOptions>());
        file.read_exact(idx_opts_slice)
            .context("Failed to read options file")?;
    }

    // Print database information
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
        mask2str(
            idx_opts.spaced_seed_mask,
            (idx_opts.l as u32)
                * if idx_opts.dna_db {
                    BITS_PER_CHAR_DNA as u32
                } else {
                    BITS_PER_CHAR_PRO as u32
                }
        )
    );
    println!("# Toggle mask = {}", mask2str(idx_opts.toggle_mask, 64));
    println!("# Total taxonomy nodes: {}", taxonomy.node_count());
    println!("# Table size: {}", kraken_index.size());
    println!("# Table capacity: {}", kraken_index.capacity());
    println!(
        "# Min clear hash value = {}",
        idx_opts.minimum_acceptable_hash_value
    );

    if idx_opts.revcom_version != CURRENT_REVCOM_VERSION as u32 {
        println!("# Built with outdated revcom version");
    }
    io::stdout().flush().context("Failed to flush stdout")?;

    // If counts are not skipped, compute and report taxonomy counts
    if !opts.skip_counts {
        let taxid_counts = kraken_index.get_value_counts();
        let mut total_seqs = 0;
        let mut taxid_counters = TaxonCountersMap::new();

        for (&taxid, &count) in taxid_counts.iter() {
            total_seqs += count;
            // Create a new TaxonCounters with the read count set to the count from the hash table
            let mut counter = TaxonCounters::new_with_precision(10);
            for _ in 0..count {
                counter.increment_read_count();
            }
            taxid_counters.insert(taxid, counter);
        }

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
                false, // report_kmer_data is always false here
                &taxonomy,
                &taxid_counters,
                total_seqs,
                0, // total_unclassified is 0
            )?;
        }
    }

    Ok(())
}
