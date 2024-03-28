use crate::{kv_store::murmurhash3, mmscanner2::MinimizerScanner};
use rayon::prelude::*;

use std::collections::HashSet;

const DEFAULT_N: usize = 4;
const RANGE_SECTIONS: usize = 1024; // must be power of 2
const RANGE_MASK: usize = RANGE_SECTIONS - 1; // must be power of 2
const MAX_N: usize = RANGE_SECTIONS;

struct EstimateCapacityOptions {
    k: usize,
    l: usize,
    n: usize,
    input_is_protein: bool,
    threads: i32,
    block_size: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
}
use clap::{Arg, Command};
use std::sync::{Arc, Mutex};

fn process_sequences(opts: &EstimateCapacityOptions) {
    let sets = Arc::new(Mutex::new(vec![HashSet::<u64>::new(); opts.n]));

    // Placeholder: Load your sequences into a collection like Vec<String>
    let sequences: Vec<String> = vec![]; // This should be populated with actual sequences

    sequences.par_iter().for_each(|seq| {
        let local_sets = sets.clone();
        process_sequence(&mut seq.clone(), opts, &mut local_sets.lock().unwrap());
    });

    // After processing, you might want to merge the sets or perform further analysis
    let final_sets = sets.lock().unwrap();
    let sum_set_sizes: usize = final_sets.iter().map(|s| s.len()).sum();
    println!("Total unique minimizers: {}", sum_set_sizes);
    // Parallel processing logic here, possibly using rayon or std::thread
}

fn process_sequence(
    seq: &mut String,
    opts: &EstimateCapacityOptions,
    sets: &mut Vec<HashSet<u64>>,
) {
    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
        0,
    );
    // Add terminator for protein sequences if not already there
    if opts.input_is_protein && seq.chars().last().unwrap() != '*' {
        seq.push('*');
    }
    scanner.load_sequence(seq, 0, 5);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.last_ambig == 1 {
            continue;
        }
        let hash_code = murmurhash3(*minimizer);
        if (hash_code & RANGE_MASK as u64) < opts.n as u64 {
            sets[(hash_code & RANGE_MASK as u64) as usize].insert(*minimizer);
        }
    }
}
/// Command line related functions
fn parse_command_line() -> EstimateCapacityOptions {
    let matches = Command::new("Kraken 2 Estimate Capacity")
        .version("1.0")
        .author("Derrick Wood <dwood@cs.jhu.edu>")
        .about("Part of the Kraken 2 taxonomic sequence classification system.")
        .arg(
            Arg::new("k")
                .long("kmer")
                .required(true)
                .help("Set length of k-mers"),
        )
        .arg(
            Arg::new("l")
                .long("minimizer")
                .required(true)
                .help("Set length of minimizers"),
        )
        .arg(
            Arg::new("n")
                .long("hashcode")
                .default_value("4")
                .help("Set maximum qualifying hash code"),
        )
        .arg(
            Arg::new("input_is_protein")
                .short('X')
                .help("Input sequences are proteins"),
        )
        .arg(
            Arg::new("threads")
                .short('p')
                .default_value("1")
                .help("Number of threads"),
        )
        .arg(
            Arg::new("block_size")
                .short('B')
                .default_value("30720000") // 30 * 1024 * 1024
                .help("Read block size"),
        )
        .arg(
            Arg::new("spaced_seed_mask")
                .short('S')
                .help("Spaced seed mask"),
        )
        .arg(
            Arg::new("toggle_mask")
                .short('T')
                .help("Minimizer ordering toggle mask"),
        )
        .get_matches();

    EstimateCapacityOptions {
        k: *matches.get_one::<usize>("k").unwrap(),
        l: *matches.get_one::<usize>("l").unwrap(),
        n: *matches.get_one::<usize>("l").unwrap(),
        input_is_protein: matches.get_flag("input_is_protein"),
        threads: *matches.get_one::<i32>("threads").unwrap(),
        block_size: *matches.get_one::<usize>("block_size").unwrap(),
        spaced_seed_mask: *matches.get_one("spaced_seed_mask").unwrap(),
        toggle_mask: *matches.get_one("toggle_mask").unwrap(),
    }
}

fn usage(exit_code: i32) {}

fn main() {
    let opts = parse_command_line(); // This would use clap or similar
    process_sequences(&opts);
}
