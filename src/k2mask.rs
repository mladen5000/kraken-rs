use rayon::prelude::*;
use std::collections::VecDeque;
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

/// Constant array to map ASCII characters to DNA bases.
/// 'A', 'C', 'G', 'T' (and their lowercase counterparts) are mapped to 0, 1, 2, 3 respectively.
/// All other characters are mapped to 4.
const ASCII2DNA: [u8; 256] = {
    let mut arr = [4; 256];
    arr[b'A' as usize] = 0;
    arr[b'C' as usize] = 1;
    arr[b'G' as usize] = 2;
    arr[b'T' as usize] = 3;
    arr[b'a' as usize] = 0;
    arr[b'c' as usize] = 1;
    arr[b'g' as usize] = 2;
    arr[b't' as usize] = 3;
    arr
};

/// Represents a DNA sequence with its header.
#[derive(Default)]
struct Sequence {
    header: String,
    seq: String,
}

/// Represents an interval of perfect matches in the sequence.
struct PerfectInterval {
    start: usize,
    finish: usize,
    left: usize,
    right: usize,
}

/// Represents a range of the sequence to be masked.
struct MaskRange {
    start: usize,
    finish: usize,
}

/// Represents the state and data required for symmetric dust masking.
struct SDust {
    kmers: VecDeque<usize>,
    perfect_intervals: Vec<PerfectInterval>,
    ranges: Vec<MaskRange>,
    seq: Sequence,
    cw: [usize; 64],
    cv: [usize; 64],
    rw: usize,
    rv: usize,
    l: usize,
}

impl SDust {
    /// Creates a new instance of `SDust` with default values.
    fn new() -> Self {
        Self {
            kmers: VecDeque::new(),
            perfect_intervals: Vec::new(),
            ranges: Vec::new(),
            seq: Sequence::default(),
            cw: [0; 64],
            cv: [0; 64],
            rw: 0,
            rv: 0,
            l: 0,
        }
    }

    /// Resets the `SDust` instance to its initial state.
    fn reset(&mut self) {
        self.kmers.clear();
        self.perfect_intervals.clear();
        self.ranges.clear();
        self.cw = [0; 64];
        self.cv = [0; 64];
        self.l = 0;
        self.rw = 0;
        self.rv = 0;
    }
}

/// Shifts the window for k-mer analysis and adjusts the counts and state accordingly.
fn shift_window(sd: &mut SDust, t: usize, window_size: usize, threshold: usize) {
    let mut s;
    if sd.kmers.len() >= window_size - 2 {
        s = sd.kmers.pop_front().unwrap();
        sd.rw -= sd.cw[s].saturating_sub(1);
        if sd.l > sd.kmers.len() {
            sd.l -= 1;
            sd.rv -= sd.cv[s].saturating_sub(1);
        }
    }
    sd.kmers.push_back(t);
    sd.l += 1;
    sd.rw += sd.cw[t];
    sd.cw[t] += 1;
    sd.rv += sd.cv[t];
    sd.cv[t] += 1;
    if sd.cv[t] * 10 > threshold * 2 {
        loop {
            s = sd.kmers[sd.kmers.len() - sd.l];
            sd.rv -= sd.cv[s].saturating_sub(1);
            sd.l -= 1;
            if s == t {
                break;
            }
        }
    }
}

/// Saves the masked regions in the sequence based on the perfect intervals detected.
fn save_masked_regions(sd: &mut SDust, window_start: usize) {
    let mut saved = false;
    if sd.perfect_intervals.is_empty() || sd.perfect_intervals.last().unwrap().start >= window_start
    {
        return;
    }
    let p = sd.perfect_intervals.last().unwrap();
    if let Some(last_range) = sd.ranges.last_mut() {
        let start = last_range.start;
        let finish = last_range.finish;
        if p.start <= finish {
            last_range.finish = std::cmp::max(p.finish, finish);
            saved = true;
        }
    }
    if !saved {
        sd.ranges.push(MaskRange {
            start: p.start,
            finish: p.finish,
        });
    }
    while !sd.perfect_intervals.is_empty()
        && sd.perfect_intervals.last().unwrap().start < window_start
    {
        sd.perfect_intervals.pop();
    }
}

/// Finds and records perfect k-mer intervals in the sequence for masking.
fn find_perfect(sd: &mut SDust, window_start: usize, threshold: usize) {
    let mut cv = sd.cv.clone();
    let mut max_left = 0;
    let mut max_right = 0;
    let mut new_left = 0;
    let mut new_right = sd.rv;

    for i in (0..sd.kmers.len() - sd.l).rev() {
        let kmer = sd.kmers[i];
        new_right += cv[kmer];
        cv[kmer] += 1;
        new_left = sd.kmers.len() - i - 1;

        if new_right * 10 > threshold * new_left {
            let mut j = 0;
            while j < sd.perfect_intervals.len()
                && sd.perfect_intervals[j].start >= i + window_start
            {
                let p = &sd.perfect_intervals[j];
                if max_right == 0 || p.right * max_left > max_right * p.left {
                    max_left = p.left;
                    max_right = p.right;
                }
                j += 1;
            }
            if max_right == 0 || new_right * max_left >= max_right * new_left {
                max_left = new_left;
                max_right = new_right;
                let p = PerfectInterval {
                    start: i + window_start,
                    finish: sd.kmers.len() + 2 + window_start,
                    left: new_left,
                    right: new_right,
                };
                sd.perfect_intervals.insert(j, p);
            }
        }
    }
}

/// Runs the symmetric dust algorithm on the sequence to detect low-complexity regions.
fn run_symmetric_dust(sd: &mut SDust, seq: &str, window_size: usize, threshold: usize) {
    let mut triplet = 0;
    let mut window_start = 0;
    let mut l = 0;

    for (i, c) in seq.chars().enumerate() {
        let base = ASCII2DNA[c as usize];
        if base < 4 {
            l += 1;
            triplet = (triplet << 2 | base as usize) & 63;
            if l >= 3 {
                window_start = std::cmp::max(l as isize - window_size as isize, 0) as usize;
                save_masked_regions(sd, window_start);
                shift_window(sd, triplet, window_size, threshold);
                if sd.rw * 10 > sd.l * threshold {
                    find_perfect(sd, window_start, threshold);
                }
            }
        }
    }
    while !sd.perfect_intervals.is_empty() {
        save_masked_regions(sd, window_start);
        window_start += 1;
    }
}

/// Prints the sequence in FASTA format to the specified output.
fn print_fasta(seq: &Sequence, out: &mut impl Write, width: usize) {
    writeln!(out, "{}", seq.header).unwrap();
    for chunk in seq.seq.as_bytes().chunks(width) {
        writeln!(out, "{}", std::str::from_utf8(chunk).unwrap()).unwrap();
    }
}

/// Processes a nucleotide to be replaced during masking.
fn process_masked_nucleotide(c: u8) -> String {
    c.to_ascii_lowercase().to_string()
}

/// The main entry point of the program.
///
/// This function processes command-line arguments, initializes data structures,
/// and runs the symmetric dust masking algorithm on input sequences.
pub fn main() {
    let args: Vec<String> = env::args().collect();
    let mut window_size = 64;
    let mut threshold = 20;
    let mut infile = String::from("/dev/stdin");
    let mut outfile = String::from("/dev/stdout");
    let mut line_width = 72;
    let mut threads = 1;

    // Start with a default closure for processing masked nucleotides
    let mut process_masked_nucleotide: Box<dyn Fn(u8) -> String + Send + Sync> =
        Box::new(|c| c.to_ascii_lowercase().to_string());

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-W" | "--window" => {
                window_size = args[i + 1].parse().unwrap();
                i += 1;
            }
            "-T" | "--level" => {
                threshold = args[i + 1].parse().unwrap();
                i += 1;
            }
            "-i" | "--in" => {
                infile = args[i + 1].clone();
                i += 1;
            }
            "-o" | "--out" => {
                outfile = args[i + 1].clone();
                i += 1;
            }
            "-w" | "--width" => {
                line_width = args[i + 1].parse().unwrap();
                i += 1;
            }
            "-t" | "--threads" => {
                threads = args[i + 1].parse().unwrap();
                i += 1;
            }
            "-r" | "--replace-masked-with" => {
                let r = args[i + 1].as_bytes()[0];
                // Assign a new closure to `process_masked_nucleotide`
                process_masked_nucleotide = Box::new(move |_| r.to_string());
                i += 1;
            }
            _ => usage(&args[0]),
        }
        i += 1;
    }

    let input = File::open(infile).unwrap();
    let mut output = File::create(outfile).unwrap();

    let sds: Arc<Mutex<Vec<_>>> =
        Arc::new(Mutex::new((0..threads).map(|_| SDust::new()).collect()));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads - 1)
        .build()
        .unwrap();
    pool.install(|| {
        let mut sds = sds.lock().unwrap();
        for sd in sds.iter_mut() {
            let sd = sd;
            mask(sd, window_size, threshold, &process_masked_nucleotide);
            print_fasta(&sd.seq, &mut output, line_width);
        }
    });
}

/// Masks the sequence using symmetric dust algorithm.
///
/// # Arguments
///
/// * `sd` - The SDust structure that holds the sequence and related data.
/// * `window_size` - The size of the window for k-mer analysis.
/// * `threshold` - The threshold for detecting low-complexity regions.
/// * `process_masked_nucleotide` - The closure for processing masked nucleotides.
fn mask<'a>(
    sd: &'a mut SDust,
    window_size: usize,
    threshold: usize,
    process_masked_nucleotide: &'a Box<dyn Fn(u8) -> String + Send + Sync>,
) -> &'a mut SDust {
    let seq_len = sd.seq.seq.len();
    let mut i = 0;

    while i < seq_len {
        if ASCII2DNA[sd.seq.seq.as_bytes()[i] as usize] != 4 {
            let start = i;

            while i < seq_len && ASCII2DNA[sd.seq.seq.as_bytes()[i] as usize] != 4 {
                let upper_case_char = sd.seq.seq[i..=i].to_ascii_uppercase(); // Clone to avoid borrowing
                sd.seq.seq.replace_range(i..=i, &upper_case_char);
                i += 1;
            }

            // Clone the sequence slice before passing to the function
            let seq_slice = sd.seq.seq[start..i].to_string(); // Clone here to avoid borrowing
            run_symmetric_dust(sd, &seq_slice, window_size, threshold);

            for range in &sd.ranges {
                for i in range.start..range.finish {
                    let processed_char = process_masked_nucleotide(sd.seq.seq[i..=i].as_bytes()[0]); // Use the closure
                    sd.seq.seq.replace_range(i..=i, &processed_char);
                }
            }

            sd.ranges.clear();
        } else {
            i += 1;
        }
    }

    sd
}

/// Displays the usage information for the program.
fn usage(prog: &str) {
    eprintln!(
        "usage: {} [-T | -level threshold] [-W | -window window size] [-i | -in input file]",
        prog
    );
    eprintln!(
        "[-w | -width sequence width] [-o | -out output file] [-r | -replace-masked-with char]"
    );
    eprintln!("[-f | -outfmt output format] [-t | -threads threads]");
    std::process::exit(1);
}
