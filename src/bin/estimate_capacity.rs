use std::collections::HashSet;
use std::io::{self, BufRead};
use std::sync::{Arc, Mutex};
use std::thread;

use kraken_rs::mmscanner::MinimizerScanner;
use kraken_rs::utilities::murmur_hash3;

const RANGE_SECTIONS: usize = 1024; // must be power of 2
const RANGE_MASK: usize = RANGE_SECTIONS - 1;
const MAX_N: usize = RANGE_SECTIONS;
const DEFAULT_N: usize = 4;
const DEFAULT_BLOCK_SIZE: usize = 30 * 1024 * 1024; // 30 MB

#[derive(Debug, Clone)]
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

impl Default for Options {
    fn default() -> Self {
        Self {
            k: 0,
            l: 0,
            n: DEFAULT_N,
            threads: 1,
            input_is_protein: false,
            spaced_seed_mask: 0,
            toggle_mask: 0,
            block_size: DEFAULT_BLOCK_SIZE,
        }
    }
}

fn parse_command_line(args: &[String]) -> Options {
    let mut opts = Options::default();
    let mut iter = args.iter().skip(1);

    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "-p" => opts.threads = iter.next().unwrap().parse().expect("Invalid thread count"),
            "-B" => opts.block_size = iter.next().unwrap().parse().expect("Invalid block size"),
            "-n" => opts.n = iter.next().unwrap().parse().expect("Invalid n value"),
            "-k" => opts.k = iter.next().unwrap().parse().expect("Invalid k value"),
            "-l" => opts.l = iter.next().unwrap().parse().expect("Invalid l value"),
            "-X" => opts.input_is_protein = true,
            "-S" => opts.spaced_seed_mask = u64::from_str_radix(iter.next().unwrap(), 2).unwrap(),
            "-T" => opts.toggle_mask = u64::from_str_radix(iter.next().unwrap(), 2).unwrap(),
            _ => usage(),
        }
    }

    if opts.k == 0 || opts.l == 0 {
        eprintln!("Missing mandatory integer parameter");
        usage();
    }
    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        usage();
    }

    opts
}

fn usage() {
    eprintln!("Usage: estimate_capacity <options>");
    std::process::exit(1);
}

fn process_sequences(opts: Options) {
    let sets: Arc<Mutex<Vec<HashSet<u64>>>> = Arc::new(Mutex::new(vec![HashSet::new(); opts.n]));
    let mut handles = vec![];

    for _ in 0..opts.threads {
        let sets_clone = Arc::clone(&sets);
        let opts_clone = opts.clone();

        let handle = thread::spawn(move || {
            let reader = io::stdin().lock().lines();
            for line in reader.flatten() {
                process_sequence(line, &opts_clone, &sets_clone);
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    let sum_set_sizes: usize = sets.lock().unwrap().iter().map(|s| s.len()).sum();
    println!("{}", (sum_set_sizes * RANGE_SECTIONS) / opts.n);
}

fn process_sequence(seq: String, opts: &Options, sets: &Arc<Mutex<Vec<HashSet<u64>>>>) {
    let mut scanner = MinimizerScanner::new(
        opts.k,
        opts.l,
        opts.spaced_seed_mask,
        !opts.input_is_protein, // rev_comp is true for DNA
        opts.toggle_mask,
        true, // enable ambiguous check
    );
    let mut seq = seq;
    if opts.input_is_protein && !seq.ends_with('*') {
        seq.push('*');
    }

    scanner.load_sequence(&seq, 0, seq.len());
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        let hash_code = murmur_hash3(minimizer);
        let index = hash_code & RANGE_MASK as u64;
        if index < opts.n as u64 {
            let mut sets_lock = sets.lock().unwrap();
            sets_lock[index as usize].insert(minimizer);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let opts = parse_command_line(&args);
    process_sequences(opts);
}
