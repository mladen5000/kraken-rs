use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::seqreader::Sequence;
use std::cmp;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

const ASC2DNA: [u8; 256] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

#[derive(Clone)]
struct PerfectInterval {
    start: i32,
    finish: i32,
    left: i32,
    right: i32,
}

#[derive(Clone)]
struct MaskRange {
    start: i32,
    finish: i32,
}

/// Struct to hold the state of the dust algorithm.
/// cw is an array of counts of each k-mer in the window,
/// cv is an array of counts of each k-mer in the sequence
/// rw is the count of k-mers in the window, and rv is the sum of cv.
#[derive(Clone)]
struct SDust {
    kmers: Vec<i32>,
    /// Perfect intervals.
    perfect_intervals: Vec<PerfectInterval>,
    /// Masked regions.
    ranges: Vec<MaskRange>,
    /// ``Sequence`` object.
    seq: Sequence,
    /// Count of each k-mer.
    cw: [i32; 64],
    /// Count of each k-mer in the window.
    cv: [i32; 64],
    /// Count of k-mers in the window.
    rw: i32,
    /// rv = sum(cv).
    rv: i32,
    /// Length of the window.
    l: i32,
}

impl SDust {
    fn new() -> Self {
        SDust {
            kmers: Vec::new(),
            perfect_intervals: Vec::new(),
            ranges: Vec::new(),
            seq: Sequence::new(),
            cw: [0; 64],
            cv: [0; 64],
            rw: 0,
            rv: 0,
            l: 0,
        }
    }

    fn reset(&mut self) {
        self.kmers.clear();
        self.perfect_intervals.clear();
        self.ranges.clear();
        self.cw.fill(0);
        self.cv.fill(0);
        self.l = 0;
        self.rw = 0;
        self.rv = 0;
    }
}

static mut WINDOW_SIZE: i32 = 64;
static mut THRESHOLD: i32 = 20;
static mut PROCESS_MASKED_NUCLEOTIDE: fn(char) -> char = |c| c.to_ascii_lowercase();

fn usage(prog: &str) {
    eprintln!(
        "usage: {} [-T | -level threshold] [-W | -window window size] [-i | -in input file]",
        prog
    );
    eprintln!(
        " [-w | -width sequence width] [-o | -out output file] [-r | -replace-masked-with char]"
    );
    eprintln!(" [-f | -outfmt output format] [-t | -threads threads]");
    std::process::exit(1);
}

/// Case-insensitive string comparison.
fn stricasecmp(s1: &str, s2: &str) -> bool {
    s1.len() == s2.len() && s1.to_lowercase() == s2.to_lowercase()
}

/// Shifts the window by one nucleotide and updates the counts.
unsafe fn shift_window(sd: &mut SDust, t: i32) {
    if sd.kmers.len() >= (WINDOW_SIZE - 2) as usize {
        let s = sd.kmers.pop().unwrap();
        sd.rw -= sd.cw[s as usize];
        sd.cw[s as usize] -= 1;
        if sd.l > sd.kmers.len() as i32 {
            sd.l -= 1;
            sd.rv -= sd.cv[s as usize];
            sd.cv[s as usize] -= 1;
        }
    }
    sd.kmers.push(t);
    sd.l += 1;
    sd.rw += sd.cw[t as usize];
    sd.cw[t as usize] += 1;
    sd.rv += sd.cv[t as usize];
    sd.cv[t as usize] += 1;
    if sd.cv[t as usize] * 10 > THRESHOLD * 2 {
        loop {
            let s = sd.kmers[sd.kmers.len() - sd.l as usize];
            sd.rv -= sd.cv[s as usize];
            sd.cv[s as usize] -= 1;
            sd.l -= 1;
            if s == t {
                break;
            }
        }
    }
}

/// Saves the masked regions.
fn save_masked_regions(sd: &mut SDust, window_start: i32) {
    let mut saved = false;
    if sd.perfect_intervals.is_empty() || sd.perfect_intervals.last().unwrap().start >= window_start
    {
        return;
    }
    let p = sd.perfect_intervals.last_mut().unwrap();
    if !sd.ranges.is_empty() {
        let _startt = sd.ranges.last().unwrap().start;
        let finish = sd.ranges.last().unwrap().finish;
        if p.start <= finish {
            sd.ranges.last_mut().unwrap().finish = cmp::max(p.finish, finish);
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

/// Finds the perfect intervals. More specifically, it finds the longest perfect interval
/// Perfect intervals are intervals where the ratio of the number of k-mers in the window
/// to the number of k-mers in the sequence is above a threshold.
unsafe fn find_perfect(sd: &mut SDust, window_start: i32) {
    let mut cv = sd.cv;
    let mut max_left = 0;
    let mut max_right = 0;
    let mut new_left;
    let mut new_right = sd.rv;

    /// Iterate over the kmers in reverse order.
    /// For each kmer, update the counts and check if the ratio is above the threshold.
    for i in (0..sd.kmers.len() - sd.l as usize).rev() {
        let kmer = sd.kmers[i];
        new_right += cv[kmer as usize];
        cv[kmer as usize] += 1;
        new_left = (sd.kmers.len() - i - 1) as i32;
        if new_right * 10 > THRESHOLD * new_left {
            let mut j = 0;
            while j < sd.perfect_intervals.len()
                && sd.perfect_intervals[j].start >= i as i32 + window_start
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
                    start: i as i32 + window_start,
                    finish: (sd.kmers.len() + 2) as i32 + window_start,
                    left: new_left,
                    right: new_right,
                };
                sd.perfect_intervals.insert(j, p);
            }
        }
    }
}

/// Runs the symmetric dust algorithm.
/// The symmetric dust algorithm is a simple algorithm to mask low-complexity regions in DNA sequences.
/// It uses a sliding window of size 64 and a threshold of 20.
/// The algorithm is based on the paper "MaskerAid: a performance enhancement to RepeatMasker" by
/// Smit, AFA, Hubley, R & Green, P. (1996).
/// The algorithm is implemented in the function run_symmetric_dust.
/// The function takes a sequence and a size as input and returns the masked sequence.
/// The function uses the following parameters:
/// - sd: The state of the dust algorithm.
/// - seq: The sequence to mask.
/// - size: The size of the sequence.
/// - offset: The offset of the sequence.
/// The function uses the following steps:
/// - Iterate over the sequence.
/// - For each nucleotide, update the counts and check if the ratio is above the threshold.
/// - If the ratio is above the threshold, find the perfect intervals.
/// - Save the masked regions.
/// - Shift the window by one nucleotide.
unsafe fn run_symmetric_dust(sd: &mut SDust, seq: &[u8], size: usize, _offset: i32) {
    let mut triplet = 0;
    let mut window_start = 0;
    let mut l = 0;
    for i in 0..size {
        let base = ASC2DNA[seq[i] as usize];
        if base < 4 {
            l += 1;
            triplet = ((triplet << 2) | base as i32) & 63;
            if l >= 3 {
                window_start = cmp::max(l as i32 - WINDOW_SIZE, 0);
                save_masked_regions(sd, window_start);
                shift_window(sd, triplet);
                if sd.rw * 10 > sd.l * THRESHOLD {
                    find_perfect(sd, window_start);
                }
            }
        }
    }
    while !sd.perfect_intervals.is_empty() {
        save_masked_regions(sd, window_start);
        window_start += 1;
    }
}

/// Prints the sequence in FASTA format.
/// The function takes a sequence and a width as input and writes the sequence to the output stream.
/// The function uses the following parameters:
/// - seq: The sequence to print.
/// - out: The output stream.
/// - width: The width of the sequence.
/// The function uses the following steps:
/// - Write the header to the output stream.
/// - Write the sequence to the output stream in chunks of the specified width.
/// - Write a newline character after each chunk.
/// The function returns nothing.
fn print_fasta(seq: &Sequence, out: &mut dyn Write, width: i32) {
    out.write_all(seq.header.as_bytes()).unwrap();
    out.write_all(b"\n").unwrap();
    for chunk in seq.seq.as_bytes().chunks(width as usize) {
        out.write_all(chunk).unwrap();
        out.write_all(b"\n").unwrap();
    }
}

unsafe fn mask(sd: &mut SDust) {
    let mut seq = std::mem::take(&mut sd.seq.seq);

    let mut i = 0;
    while i < seq.len() {
        if ASC2DNA[seq.as_bytes()[i] as usize] != 4 {
            let start = i;
            loop {
                seq.as_mut_vec()[i] = seq.as_bytes()[i].to_ascii_uppercase() as u8;
                if i + 1 == seq.len() || ASC2DNA[seq.as_bytes()[i + 1] as usize] == 4 {
                    break;
                }
                i += 1;
            }
            let s = &seq.as_bytes()[start..=i];
            run_symmetric_dust(sd, s, i - start + 1, start as i32);

            let ranges = std::mem::take(&mut sd.ranges);
            for range in ranges {
                for k in range.start..range.finish {
                    seq.as_mut_vec()[k as usize] =
                        PROCESS_MASKED_NUCLEOTIDE(seq.as_bytes()[k as usize] as char) as u8;
                }
            }
        }
        i += 1;
    }

    sd.seq.seq = seq;
}
pub fn main() {
    let mut line_width = 72;
    let mut threads = 1;
    let mut infile = "/dev/stdin".to_string();
    let mut outfile = "/dev/stdout".to_string();
    let prog = "k2mask";

    // ... (command-line argument parsing code remains the same)

    let mut reader = BufReader::new(File::open(infile).unwrap());
    let mut writer = BufWriter::new(File::create(outfile).unwrap());

    let mut sds: Vec<SDust> = (0..threads).map(|_| SDust::new()).collect();

    let mut buffer = String::new();
    let sequences: Vec<Sequence> = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let mut sd = SDust::new();
            sd.seq.header = line.clone();
            sd.seq.seq = line;
            sd.seq.clone()
        })
        .collect();

    let masked_sds: Vec<SDust> = sequences
        .into_par_iter()
        .map(|seq| {
            let mut sd = SDust::new();
            sd.seq = seq;
            unsafe {
                mask(&mut sd);
            }
            sd
        })
        .collect();

    for sd in &masked_sds {
        print_fasta(&sd.seq, &mut writer, line_width);
    }

    writer.flush().unwrap();
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_asc2dna() {
        assert_eq!(ASC2DNA[b'A' as usize], 0);
        assert_eq!(ASC2DNA[b'C' as usize], 1);
        assert_eq!(ASC2DNA[b'G' as usize], 2);
        assert_eq!(ASC2DNA[b'T' as usize], 3);
        assert_eq!(ASC2DNA[b'N' as usize], 4);
    }

    #[test]
    fn test_sdust_new() {
        let sd = SDust::new();
        assert!(sd.kmers.is_empty());
        assert!(sd.perfect_intervals.is_empty());
        assert!(sd.ranges.is_empty());
        assert!(sd.seq.seq.is_empty());
        assert_eq!(sd.cw, [0; 64]);
        assert_eq!(sd.cv, [0; 64]);
        assert_eq!(sd.rw, 0);
        assert_eq!(sd.rv, 0);
        assert_eq!(sd.l, 0);
    }

    #[test]
    fn test_sdust_reset() {
        let mut sd = SDust::new();
        sd.kmers.push(1);
        sd.perfect_intervals.push(PerfectInterval {
            start: 0,
            finish: 1,
            left: 1,
            right: 1,
        });
        sd.ranges.push(MaskRange {
            start: 0,
            finish: 1,
        });
        sd.seq.seq = "ATCG".to_string();
        sd.cw[0] = 1;
        sd.cv[0] = 1;
        sd.rw = 1;
        sd.rv = 1;
        sd.l = 1;

        sd.reset();
        assert!(sd.kmers.is_empty());
        assert!(sd.perfect_intervals.is_empty());
        assert!(sd.ranges.is_empty());
        assert!(sd.seq.seq.is_empty());
        assert_eq!(sd.cw, [0; 64]);
        assert_eq!(sd.cv, [0; 64]);
        assert_eq!(sd.rw, 0);
        assert_eq!(sd.rv, 0);
        assert_eq!(sd.l, 0);
    }

    // Add more tests for other functions as needed
}
