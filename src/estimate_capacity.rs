/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::HashSet;
use std::io::{self, Read};
use std::sync::{Arc, Mutex};

use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;

const RANGE_SECTIONS: usize = 1024; // Must be a power of 2
const RANGE_MASK: usize = RANGE_SECTIONS - 1;
const MAX_N: usize = RANGE_SECTIONS;
const DEFAULT_N: usize = 4;
const DEFAULT_BLOCK_SIZE: usize = 30 * 1024 * 1024; // 30 MB

const DEFAULT_SPACED_SEED_MASK: u64 = 0xFFFFFFFFFFFFFFFF;
const DEFAULT_TOGGLE_MASK: u64 = 0x0;

#[derive(Clone)]
struct Options {
    k: usize,
    l: usize,
    n: usize,
    input_is_protein: bool,
    threads: usize,
    block_size: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
}

#[derive(Clone)]
struct Sequence {
    seq: String,
}

struct BatchSequenceReader<R: Read> {
    reader: R,
    buffer: Vec<u8>,
    position: usize,
    capacity: usize,
}

impl<R: Read> BatchSequenceReader<R> {
    fn new(reader: R, capacity: usize) -> Self {
        Self {
            reader,
            buffer: Vec::with_capacity(capacity),
            position: 0,
            capacity,
        }
    }

    fn load_block(&mut self) -> io::Result<bool> {
        self.buffer.clear();
        self.buffer.resize(self.capacity, 0);
        let bytes_read = self.reader.read(&mut self.buffer)?;
        self.buffer.truncate(bytes_read);
        self.position = 0;
        Ok(bytes_read > 0)
    }

    fn next_sequence(&mut self) -> Option<Sequence> {
        if self.position >= self.buffer.len() {
            return None;
        }
        let start = self.position;
        while self.position < self.buffer.len() && self.buffer[self.position] != b'\n' {
            self.position += 1;
        }
        let end = self.position;
        if self.position < self.buffer.len() {
            self.position += 1; // Skip the newline
        }
        let seq = String::from_utf8_lossy(&self.buffer[start..end]).to_string();
        Some(Sequence { seq })
    }
}

fn parse_command_line() -> Options {
    let matches = Command::new("estimate_capacity")
        .version("1.0")
        .about("Estimate the capacity needed for Kraken 2 database")
        .arg(
            Arg::new("k")
                .short('k')
                .long("k")
                .required(true)
                .help("Set length of k-mers")
                .value_name("K_VALUE")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("l")
                .short('l')
                .long("l")
                .required(true)
                .help("Set length of minimizers")
                .value_name("L_VALUE")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("n")
                .short('n')
                .long("n")
                .help("Set maximum qualifying hash code")
                .value_name("N_VALUE")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("input_is_protein")
                .short('X')
                .long("input_is_protein")
                .help("Input sequences are proteins")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("spaced_seed_mask")
                .short('S')
                .long("spaced_seed_mask")
                .help("Spaced seed mask (binary string)")
                .value_name("SPACED_SEED_MASK")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("toggle_mask")
                .short('T')
                .long("toggle_mask")
                .help("Minimizer ordering toggle mask (binary string)")
                .value_name("TOGGLE_MASK")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("block_size")
                .short('B')
                .long("block_size")
                .help("Read block size")
                .value_name("BLOCK_SIZE")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("threads")
                .short('p')
                .long("threads")
                .help("Number of threads")
                .value_name("THREADS")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .get_matches();

    // Access and parse arguments
    let k = *matches
        .get_one::<usize>("k")
        .expect("k is required and must be a positive integer");

    let l = *matches
        .get_one::<usize>("l")
        .expect("l is required and must be a positive integer");

    let n = *matches.get_one::<usize>("n").unwrap_or(&DEFAULT_N);

    if n > MAX_N {
        eprintln!("Error: n must be no more than {}", MAX_N);
        std::process::exit(1);
    }

    let input_is_protein = matches.get_flag("input_is_protein");

    let spaced_seed_mask = matches
        .get_one::<String>("spaced_seed_mask")
        .map(|v| {
            u64::from_str_radix(v, 2).unwrap_or_else(|_| {
                eprintln!("Error: spaced_seed_mask must be a binary string");
                std::process::exit(1);
            })
        })
        .unwrap_or(DEFAULT_SPACED_SEED_MASK);

    let toggle_mask = matches
        .get_one::<String>("toggle_mask")
        .map(|v| {
            u64::from_str_radix(v, 2).unwrap_or_else(|_| {
                eprintln!("Error: toggle_mask must be a binary string");
                std::process::exit(1);
            })
        })
        .unwrap_or(DEFAULT_TOGGLE_MASK);

    let block_size = *matches
        .get_one::<usize>("block_size")
        .unwrap_or(&DEFAULT_BLOCK_SIZE);

    let threads = *matches.get_one::<usize>("threads").unwrap_or(&1);

    if k < l {
        eprintln!("Error: k cannot be less than l");
        std::process::exit(1);
    }

    Options {
        k,
        l,
        n,
        input_is_protein,
        threads,
        block_size,
        spaced_seed_mask,
        toggle_mask,
    }
}
fn main() {
    let opts = parse_command_line();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build_global()
        .unwrap();

    process_sequences(&opts);
}

fn process_sequences(opts: &Options) {
    let sets: Arc<Vec<Mutex<HashSet<u64>>>> =
        Arc::new((0..opts.n).map(|_| Mutex::new(HashSet::new())).collect());

    let stdin = io::stdin();
    let mut reader = BatchSequenceReader::new(stdin.lock(), opts.block_size);

    while reader.load_block().unwrap() {
        let sequences: Vec<Sequence> = std::iter::from_fn(|| reader.next_sequence()).collect();

        sequences.par_iter().for_each(|sequence| {
            process_sequence(&sequence.seq, opts, &sets);
        });
    }

    let sum_set_sizes: usize = sets
        .iter()
        .map(|set_mutex| set_mutex.lock().unwrap().len())
        .sum();

    let sum_set_sizes = sum_set_sizes + 1; // Ensure non-zero estimate

    println!("{}", (sum_set_sizes * RANGE_SECTIONS) / opts.n);
}

fn process_sequence(seq: &str, opts: &Options, sets: &Arc<Vec<Mutex<HashSet<u64>>>>) {
    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
        1, // Assuming revcom_version is 1
    );

    let mut seq = seq.to_string();

    // Add terminator for protein sequences if not already there
    if opts.input_is_protein && !seq.ends_with('*') {
        seq.push('*');
    }

    scanner.load_sequence(&seq);

    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        let hash_code = murmur_hash3(minimizer);

        if (hash_code & RANGE_MASK as u64) < opts.n as u64 {
            let idx = (hash_code & RANGE_MASK as u64) as usize;
            let mut set = sets[idx].lock().unwrap();
            set.insert(minimizer);
        }
    }
}

struct MinimizerScanner {
    sequence: Vec<char>,
    position: usize,
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    is_dna: bool,
    toggle_mask: u64,
}

impl MinimizerScanner {
    fn new(k: usize, l: usize, spaced_seed_mask: u64, is_dna: bool, toggle_mask: u64) -> Self {
        Self {
            sequence: Vec::new(),
            position: 0,
            k,
            l,
            spaced_seed_mask,
            is_dna,
            toggle_mask,
        }
    }

    fn load_sequence(&mut self, seq: &str) {
        self.sequence = seq.chars().collect();
        self.position = 0;
    }

    fn next_minimizer(&mut self) -> Option<u64> {
        // Placeholder implementation: returns k-mers as minimizers
        if self.position + self.k <= self.sequence.len() {
            let kmer: String = self.sequence[self.position..self.position + self.k]
                .iter()
                .collect();
            self.position += 1;
            let minimizer = kmer_to_u64(&kmer);
            Some(minimizer)
        } else {
            None
        }
    }

    fn is_ambiguous(&self) -> bool {
        // Placeholder: Always returns false
        false
    }
}

fn kmer_to_u64(kmer: &str) -> u64 {
    // Convert k-mer string to u64
    let mut result = 0u64;
    for c in kmer.chars() {
        let val = match c {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => 0, // Treat other characters as 0
        };
        result = (result << 2) | val;
    }
    result
}

fn murmur_hash3(key: u64) -> u64 {
    // Placeholder for MurmurHash3 function
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write_u64(key);
    hasher.finish()
}
