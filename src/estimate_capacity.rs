use std::collections::{HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;

const RANGE_SECTIONS: usize = 1024; // must be power of 2
const RANGE_MASK: usize = RANGE_SECTIONS - 1;
const MAX_N: usize = RANGE_SECTIONS;
const DEFAULT_N: usize = 4;
const DEFAULT_BLOCK_SIZE: usize = 30 * 1024 * 1024; // yes, 30 MB

struct Options {
    k: usize,
    l: usize,
    n: usize,
    input_is_protein: bool,
    threads: i32,
    block_size: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
}

fn parse_command_line(args: &[String]) -> Options {
    let mut opts = Options {
        k: 0,
        l: 0,
        n: DEFAULT_N,
        input_is_protein: false,
        threads: 1,
        block_size: DEFAULT_BLOCK_SIZE,
        spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
        toggle_mask: DEFAULT_TOGGLE_MASK,
    };

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => usage(),
            "-p" => {
                let sig = parse_arg_as_i32(&args, i);
                if sig < 1 {
                    eprintln!("must have at least 1 thread");
                    usage();
                }
                opts.threads = sig;
                i += 2;
            }
            "-B" => {
                let sig = parse_arg_as_usize(&args, i);
                if sig < 1 {
                    eprintln!("block size must be positive");
                    usage();
                }
                opts.block_size = sig;
                i += 2;
            }
            "-n" => {
                let sig = parse_arg_as_usize(&args, i);
                if sig < 1 {
                    eprintln!("n must be positive integer");
                    usage();
                }
                if sig > MAX_N {
                    eprintln!("n must be no more than {}", MAX_N);
                    usage();
                }
                opts.n = sig;
                i += 2;
            }
            "-k" => {
                let sig = parse_arg_as_usize(&args, i);
                if sig < 1 {
                    eprintln!("k must be positive integer");
                    usage();
                }
                opts.k = sig;
                i += 2;
            }
            "-l" => {
                let sig = parse_arg_as_usize(&args, i);
                if sig < 1 {
                    eprintln!("l must be positive integer");
                    usage();
                }
                if sig > 31 {
                    eprintln!("l must be no more than 31");
                    usage();
                }
                opts.l = sig;
                i += 2;
            }
            "-X" => {
                opts.input_is_protein = true;
                i += 1;
            }
            "-S" => {
                opts.spaced_seed_mask = parse_arg_as_u64(&args, i);
                i += 2;
            }
            "-T" => {
                opts.toggle_mask = parse_arg_as_u64(&args, i);
                i += 2;
            }
            _ => {
                eprintln!("Invalid option: {}", args[i]);
                usage();
            }
        }
    }

    if opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        expand_spaced_seed_mask(
            &mut opts.spaced_seed_mask,
            if opts.input_is_protein {
                BITS_PER_CHAR_PRO
            } else {
                BITS_PER_CHAR_DNA
            },
        );
    }

    if opts.k == 0 || opts.l == 0 {
        eprintln!("missing mandatory integer parameter");
        usage();
    }

    opts
}

fn parse_arg_as_usize(args: &[String], index: usize) -> usize {
    args[index + 1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid argument: {}", args[index + 1]);
        usage();
    })
}

fn parse_arg_as_i32(args: &[String], index: usize) -> i32 {
    args[index + 1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid argument: {}", args[index + 1]);
        usage();
    })
}

fn parse_arg_as_u64(args: &[String], index: usize) -> u64 {
    args[index + 1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid argument: {}", args[index + 1]);
        usage();
    })
}

fn usage() {
    eprintln!(
        "Usage: estimate_capacity <options>\n\
    \n\
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
    process::exit(1);
}

fn expand_spaced_seed_mask(spaced_seed_mask: &mut u64, bits_per_char: usize) {
    let mask = (1 << bits_per_char) - 1;
    let mut new_mask = 0;
    for i in 0..bits_per_char {
        if (*spaced_seed_mask & (1 << i)) != 0 {
            new_mask |= mask << i;
        }
    }
    *spaced_seed_mask = new_mask;
}

fn process_sequence(seq: &str, opts: &Options, sets: &mut Vec<HashSet<u64>>) {
    let mut scanner = MinimizerScanner::new(
        opts.k,
        opts.l,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
    );
    let seq = if opts.input_is_protein {
        format!("{}*", seq)
    } else {
        seq.to_string()
    };
    scanner.load_sequence(&seq);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        let hash_code = murmur_hash3(minimizer);
        if (hash_code & RANGE_MASK) < opts.n {
            let index = (hash_code & RANGE_MASK) as usize;
            sets[index].insert(minimizer);
        }
    }
}

fn process_sequences(opts: &Options) -> usize {
    let mut sets: Vec<HashSet<u64>> = vec![HashSet::new(); opts.n];

    let stdin = io::stdin();
    let reader = BufReader::new(stdin.lock());

    let mut sum_set_sizes = 0;
    let mut have_work = true;
    let mut reader = BatchSequenceReader::new();
    let mut sequence = Sequence::new();

    while have_work {
        have_work = reader.load_block(&mut reader, opts.block_size);
        if have_work {
            while reader.next_sequence(&mut sequence) {
                process_sequence(&sequence.seq, opts, &mut sets);
            }
        }
    }

    for s in &sets {
        sum_set_sizes += s.len();
    }
    sum_set_sizes += 1; // ensure non-zero estimate
    (sum_set_sizes as f64 * RANGE_SECTIONS as f64 / opts.n as f64) as usize
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let opts = parse_command_line(&args);
    process::exit(process_sequences(&opts));
}
