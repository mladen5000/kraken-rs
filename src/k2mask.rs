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

#[derive(Clone)]
struct SDust {
    kmers: Vec<i32>,
    perfect_intervals: Vec<PerfectInterval>,
    ranges: Vec<MaskRange>,
    seq: Sequence,
    cw: [i32; 64],
    cv: [i32; 64],
    rw: i32,
    rv: i32,
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

fn stricasecmp(s1: &str, s2: &str) -> bool {
    s1.len() == s2.len() && s1.to_lowercase() == s2.to_lowercase()
}

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

fn save_masked_regions(sd: &mut SDust, window_start: i32) {
    let mut saved = false;
    if sd.perfect_intervals.is_empty() || sd.perfect_intervals.last().unwrap().start >= window_start
    {
        return;
    }
    let p = sd.perfect_intervals.last_mut().unwrap();
    if !sd.ranges.is_empty() {
        let start = sd.ranges.last().unwrap().start;
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

unsafe fn find_perfect(sd: &mut SDust, window_start: i32) {
    let mut cv = sd.cv;
    let mut max_left = 0;
    let mut max_right = 0;
    let mut new_left;
    let mut new_right = sd.rv;

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

unsafe fn run_symmetric_dust(sd: &mut SDust, seq: &[u8], size: usize, offset: i32) {
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

fn print_fasta(seq: &Sequence, out: &mut dyn Write, width: i32) {
    out.write_all(seq.header.as_bytes()).unwrap();
    out.write_all(b"\n").unwrap();
    for chunk in seq.seq.as_bytes().chunks(width as usize) {
        out.write_all(chunk).unwrap();
        out.write_all(b"\n").unwrap();
    }
}

fn mask(sd: &mut SDust) {
    let mut i = 0;
    let seq = &mut sd.seq.seq;

    while i < seq.len() {
        if ASC2DNA[seq.as_bytes()[i] as usize] != 4 {
            let start = i;
            loop {
                unsafe { seq.as_bytes_mut()[i] = seq.as_bytes()[i].to_ascii_uppercase() }
                if i + 1 == seq.len() || ASC2DNA[seq.as_bytes()[i + 1] as usize] == 4 {
                    break;
                }
                i += 1;
            }
            let s = unsafe { &mut seq.as_bytes_mut()[start..=i] };
            unsafe { run_symmetric_dust(sd, s, i - start + 1, start as i32) };
            for j in 0..sd.ranges.len() {
                for k in sd.ranges[j].start..sd.ranges[j].finish {
                    unsafe {
                        s[k as usize] = PROCESS_MASKED_NUCLEOTIDE(s[k as usize] as char) as u8
                    };
                }
            }
            sd.ranges.clear();
        }
        i += 1;
    }
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
