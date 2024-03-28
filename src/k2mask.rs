use std::cmp;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::sync::{Arc, Mutex};

use flate2::read::GzDecoder;
use threadpool::ThreadPool;

use crate::seqreader::*;

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

#[derive(Default, Debug, Clone)]
struct PerfectInterval {
    start: usize,
    finish: usize,
    left: usize,
    right: usize,
}

#[derive(Default, Debug, Clone)]
struct MaskRange {
    start: usize,
    finish: usize,
}

#[derive(Debug, Clone)]
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
    fn new() -> Self {
        SDust {
            kmers: VecDeque::new(),
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
        self.cw = [0; 64];
        self.cv = [0; 64];
        self.l = 0;
        self.rw = 0;
        self.rv = 0;
    }
}

fn stricasecmp(s1: &str, s2: &str) -> bool {
    if s1.len() != s2.len() {
        return false;
    }
    s1.chars()
        .zip(s2.chars())
        .all(|(c1, c2)| c1.to_ascii_lowercase() == c2.to_ascii_lowercase())
}

fn shift_window(sd: &mut SDust, t: usize, window_size: usize) {
    if sd.kmers.len() >= window_size - 2 {
        let s = sd.kmers.pop_front().unwrap();
        sd.rw -= sd.cw[s];
        sd.cw[s] -= 1;
        if sd.l > sd.kmers.len() {
            sd.l -= 1;
            sd.rv -= sd.cv[s];
            sd.cv[s] -= 1;
        }
    }
    sd.kmers.push_back(t);
    sd.l += 1;
    sd.rw += sd.cw[t];
    sd.cw[t] += 1;
    sd.rv += sd.cv[t];
    sd.cv[t] += 1;
    if sd.cv[t] * 10 > 20 * 2 {
        loop {
            let s = sd.kmers[sd.kmers.len() - sd.l];
            sd.rv -= sd.cv[s];
            sd.cv[s] -= 1;
            sd.l -= 1;
            if s == t {
                break;
            }
        }
    }
}

fn save_masked_regions(sd: &mut SDust, window_start: usize) {
    let mut saved = false;
    if sd.perfect_intervals.is_empty() || sd.perfect_intervals.last().unwrap().start >= window_start
    {
        return;
    }
    let p = sd.perfect_intervals.last().unwrap();
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

fn find_perfect(sd: &mut SDust, window_start: usize, threshold: usize) {
    let mut cv = [0; 64];
    let mut max_left = 0;
    let mut max_right = 0;
    let mut new_left;
    let mut new_right = sd.rv;

    cv.copy_from_slice(&sd.cv);
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

fn run_symmetric_dust(
    sd: &mut SDust,
    seq: &[u8],
    offset: usize,
    window_size: usize,
    threshold: usize,
) {
    let mut triplet = 0;
    let mut window_start = 0;
    let mut l = 0;
    for &base in seq {
        let base = ASC2DNA[base as usize];
        if base < 4 {
            l += 1;
            triplet = ((triplet << 2) | base as usize) & 63;
            if l >= 3 {
                window_start = cmp::max(l as isize - window_size as isize, 0) as usize;
                save_masked_regions(sd, window_start);
                shift_window(sd, triplet, window_size);
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

fn print_fasta<W: Write>(seq: &Sequence, out: &mut W, width: usize) -> io::Result<()> {
    writeln!(out, "{}", seq.header)?;
    for chunk in seq.seq.as_bytes().chunks(width) {
        writeln!(out, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

fn mask(sd: &mut SDust, process_masked_nucleotide: &dyn Fn(char) -> char) {
    let seq = sd.seq.seq.as_bytes();
    let mut i = 0;
    while i < seq.len() {
        if ASC2DNA[seq[i] as usize] != 4 {
            let start = i;
            while i + 1 < seq.len() && ASC2DNA[seq[i + 1] as usize] != 4 {
                i += 1;
            }
            let s = unsafe { &mut sd.seq.seq.as_bytes_mut()[start..=i] };
            run_symmetric_dust(sd, s, start, 64, 20);
            for range in &sd.ranges {
                for j in range.start..range.finish {
                    s[j - start] = process_masked_nucleotide(s[j - start] as char) as u8;
                }
            }
            sd.ranges.clear();
        }
        i += 1;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let prog = &args[0];

    let mut window_size = 64;
    let mut threshold = 20;
    let mut infile = "/dev/stdin".to_string();
    let mut outfile = "/dev/stdout".to_string();
    let mut line_width = 72;
    let mut threads = 1;
    let mut process_masked_nucleotide: Box<dyn Fn(char) -> char> =
        Box::new(|c| c.to_ascii_lowercase());

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-W" | "--window" => {
                window_size = args[i + 1].parse().expect("Invalid window size");
                i += 1;
            }
            "-T" | "--level" => {
                threshold = args[i + 1].parse().expect("Invalid threshold");
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
                line_width = args[i + 1].parse().expect("Invalid line width");
                i += 1;
            }
            "-f" | "--outfmt" => {
                let format = args[i + 1].clone();
                if !stricasecmp(&format, "fasta") {
                    eprintln!("{}: currently only supports outputting FASTA.", prog);
                    std::process::exit(1);
                }
                i += 1;
            }
            "-t" | "--threads" => {
                threads = args[i + 1].parse().expect("Invalid number of threads");
                i += 1;
            }
            "-r" | "--replace-masked-with" => {
                let r = args[i + 1]
                    .chars()
                    .next()
                    .expect("Invalid replacement character");
                process_masked_nucleotide = Box::new(move |_| r);
                i += 1;
            }
            "-h" | "--help" => {
                eprintln!("usage: {} [-T | -level threshold] [-W | -window window size] [-i | -in input file]", prog);
                eprintln!("       [-w | -width sequence width] [-o | -out output file] [-r | -replace-masked-with char]");
                eprintln!("       [-f | -outfmt output format] [-t | -threads threads]");
                std::process::exit(1);
            }
            _ => {
                eprintln!("{}: reads from stdin and writes to stdout", prog);
                eprintln!("usage: {} [-T | -level threshold] [-W | -window window size] [-i | -in input file]", prog);
                eprintln!("       [-w | -width sequence width] [-o | -out output file] [-r | -replace-masked-with char]");
                eprintln!("       [-f | -outfmt output format] [-t | -threads threads]");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let infile = File::open(infile).expect("Failed to open input file");
    let infile = BufReader::new(GzDecoder::new(infile));
    let outfile = File::create(outfile).expect("Failed to create output file");
    let mut outfile = io::BufWriter::new(outfile);

    let pool = ThreadPool::new(threads);
    let sds: Arc<Mutex<Vec<SDust>>> = Arc::new(Mutex::new(vec![SDust::new(); threads]));

    let mut reader = BatchSequenceReader::new();
    while let Some(seq) = reader.read_next_sequence() {
        let mut sd = sds.lock().unwrap().pop().unwrap();
        sd.seq = seq;
        let sds_clone = Arc::clone(&sds);
        let process_masked_nucleotide = process_masked_nucleotide.clone();
        pool.execute(move || {
            mask(&mut sd, &process_masked_nucleotide);
            sds_clone.lock().unwrap().push(sd);
        });
    }

    pool.join();

    let sds = Arc::try_unwrap(sds).unwrap().into_inner().unwrap();
    for sd in sds {
        print_fasta(&sd.seq, &mut outfile, line_width).expect("Failed to write output");
    }
}
