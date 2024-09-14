use rayon::prelude::*;
use std::collections::HashSet;
use std::env;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

const RANGE_SECTIONS: usize = 1024;
const RANGE_MASK: usize = RANGE_SECTIONS - 1;
const MAX_N: usize = RANGE_SECTIONS;
const DEFAULT_N: usize = 4;
const DEFAULT_BLOCK_SIZE: usize = 30 * 1024 * 1024;
const DEFAULT_SPACED_SEED_MASK: u64 = 0;
const DEFAULT_TOGGLE_MASK: u64 = 0;

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

fn main() {
    let mut opts = Options {
        k: 0,
        l: 0,
        n: DEFAULT_N,
        threads: 1,
        input_is_protein: false,
        spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
        toggle_mask: DEFAULT_TOGGLE_MASK,
        block_size: DEFAULT_BLOCK_SIZE,
    };
    parse_command_line(env::args().collect(), &mut opts);
    process_sequences(&opts);
}

fn parse_command_line(args: Vec<String>, opts: &mut Options) {
    let mut args_iter = args.iter();
    args_iter.next(); // Skip the program name

    while let Some(arg) = args_iter.next() {
        match arg.as_str() {
            "-h" | "?" => usage(0),
            "-p" => {
                opts.threads = args_iter
                    .next()
                    .expect("Expected a value for threads")
                    .parse()
                    .expect("Invalid thread count");
                if opts.threads < 1 {
                    eprintln!("must have at least 1 thread");
                    usage(1);
                }
            }
            "-B" => {
                opts.block_size = args_iter
                    .next()
                    .expect("Expected a value for block size")
                    .parse()
                    .expect("Invalid block size");
                if opts.block_size < 1 {
                    eprintln!("block size must be positive");
                    usage(1);
                }
            }
            "-n" => {
                opts.n = args_iter
                    .next()
                    .expect("Expected a value for n")
                    .parse()
                    .expect("Invalid value for n");
                if opts.n < 1 {
                    eprintln!("n must be positive integer");
                    usage(1);
                }
                if opts.n > MAX_N {
                    eprintln!("n must be no more than {}", MAX_N);
                    usage(1);
                }
            }
            "-k" => {
                opts.k = args_iter
                    .next()
                    .expect("Expected a value for k")
                    .parse()
                    .expect("Invalid value for k");
                if opts.k < 1 {
                    eprintln!("k must be positive integer");
                    usage(1);
                }
            }
            "-l" => {
                opts.l = args_iter
                    .next()
                    .expect("Expected a value for l")
                    .parse()
                    .expect("Invalid value for l");
                if opts.l < 1 {
                    eprintln!("l must be positive integer");
                    usage(1);
                }
                if opts.l > 31 {
                    eprintln!("l must be no more than 31");
                    usage(1);
                }
            }
            "-X" => opts.input_is_protein = true,
            "-S" => {
                opts.spaced_seed_mask = u64::from_str_radix(
                    args_iter
                        .next()
                        .expect("Expected a value for spaced seed mask"),
                    2,
                )
                .expect("Invalid spaced seed mask");
            }
            "-T" => {
                opts.toggle_mask = u64::from_str_radix(
                    args_iter.next().expect("Expected a value for toggle mask"),
                    2,
                )
                .expect("Invalid toggle mask");
            }
            _ => {
                eprintln!("Unknown option: {}", arg);
                usage(1);
            }
        }
    }

    if opts.k == 0 || opts.l == 0 {
        eprintln!("missing mandatory integer parameter");
        usage(1);
    }

    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        usage(1);
    }
}

fn process_sequences(opts: &Options) {
    let sets = Arc::new(Mutex::new(vec![HashSet::<u64>::new(); opts.n]));

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build_global()
        .unwrap();

    (0..opts.threads).into_par_iter().for_each(|_| {
        let sets = Arc::clone(&sets);
        let mut reader = BatchSequenceReader::new();
        let mut sequence = String::new();
        let mut have_work = true;

        while have_work {
            {
                let mut sets = sets.lock().unwrap();
                have_work = reader.load_block(std::io::stdin().lock(), opts.block_size);
                if have_work {
                    while reader.next_sequence(&mut sequence) {
                        process_sequence(&mut sequence, opts, &mut sets);
                    }
                }
            }
            thread::sleep(Duration::from_millis(10)); // allow other threads to acquire lock
        }
    });

    let sum_set_sizes: usize = sets
        .lock()
        .unwrap()
        .iter()
        .map(|set| set.len())
        .sum::<usize>()
        + 1; // ensure non-zero estimate

    println!(
        "{}",
        ((sum_set_sizes as f64) * RANGE_SECTIONS as f64 * 1.0 / opts.n as f64) as usize
    );
}

fn process_sequence(seq: &mut String, opts: &Options, sets: &mut Vec<HashSet<u64>>) {
    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
    );

    if opts.input_is_protein && !seq.ends_with('*') {
        seq.push('*');
    }

    scanner.load_sequence(seq);

    while let Some(minimizer) = scanner.next_minimizer().cloned() {
        let is_ambiguous = scanner.is_ambiguous(); // Now it's safe to use the immutable borrow

        if is_ambiguous {
            continue;
        }

        let hash_code = murmurhash3(minimizer); // Use the cloned minimizer
        if (hash_code & RANGE_MASK as u64) < opts.n as u64 {
            sets[(hash_code & RANGE_MASK as u64) as usize].insert(minimizer);
        }
    }
}

fn usage(exit_code: i32) {
    eprintln!(
        "Usage: estimate_capacity <options>\n\n\
         Options (*mandatory):\n\
         * -k INT        Set length of k-mers\n\
         * -l INT        Set length of minimizers\n\
         -n INT        Set maximum qualifying hash code\n\
         -X            Input sequences are proteins\n\
         -S BITSTRING  Spaced seed mask\n\
         -T BITSTRING  Minimizer ordering toggle mask\n\
         -B INT        Read block size\n\
         -p INT        Number of threads"
    );
    std::process::exit(exit_code);
}

// You'll need to implement these functions or import them from other modules
fn murmurhash3(key: u64) -> u64 {
    // Implement MurmurHash3 here
    unimplemented!()
}

struct BatchSequenceReader;
struct MinimizerScanner;

// Implement BatchSequenceReader and MinimizerScanner as needed
