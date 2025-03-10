use std::string::String;
use std::io;

const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
const DEFAULT_SPACED_SEED_MASK: u64 = 0;
const BITS_PER_CHAR_DNA: usize = 2;
const BITS_PER_CHAR_PRO: usize = 4;
const CURRENT_REVCOM_VERSION: i32 = 1;

#[derive(Clone)]
struct MinimizerData {
    candidate: u64,
    pos: isize,
}

pub struct MinimizerScanner {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    rev_comp: bool,
    toggle_mask: u64,
    ambiguous_check: bool,
    current_minimizer: u64,
    is_ambiguous: bool,
    sequence: String,
    pos: usize,
}

impl MinimizerScanner {
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        rev_comp: bool,
        toggle_mask: u64,
        ambiguous_check: bool,
    ) -> Self {
        Self {
            k,
            l,
            spaced_seed_mask,
            rev_comp,
            toggle_mask,
            ambiguous_check,
            current_minimizer: 0,
            is_ambiguous: false,
            sequence: String::new(),
            pos: 0,
        }
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, end: usize) {
        self.sequence = seq[start..end].to_string();
        self.pos = 0;
        self.is_ambiguous = false;
        self.current_minimizer = 0;
    }
    
    pub fn load_sequence_with_range(&mut self, seq: &str, start: usize, end: usize) {
        // This is an alias to load_sequence for compatibility with the build_db.rs calls
        self.load_sequence(seq, start, end);
    }

    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.pos >= self.sequence.len() {
            return None;
        }

        let kmer_end = self.pos + self.k;
        if kmer_end > self.sequence.len() {
            self.pos = self.sequence.len();
            return None;
        }

        let kmer = &self.sequence[self.pos..kmer_end];
        self.pos += 1;

        if self.ambiguous_check {
            self.is_ambiguous = kmer.chars().any(|c| match c {
                'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't' => false,
                _ => true
            });
            if self.is_ambiguous {
                return Some(0);
            }
        }

        let mut value = self.calculate_minimizer(kmer);
        if self.rev_comp {
            let rc_value = self.calculate_minimizer(&revcomp(kmer));
            value = std::cmp::min(value, rc_value);
        }
        
        self.current_minimizer = value;
        Some(value)
    }

    pub fn is_ambiguous(&self) -> bool {
        self.is_ambiguous
    }

    pub fn last_minimizer(&self) -> u64 {
        self.current_minimizer
    }

    fn calculate_minimizer(&self, kmer: &str) -> u64 {
        let mut value = 0u64;
        for (i, c) in kmer.chars().enumerate() {
            if (self.spaced_seed_mask & (1 << i)) != 0 {
                value = (value << 2) | match c.to_ascii_uppercase() {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    _ => 0,
                };
            }
        }
        if (self.toggle_mask & value) != 0 {
            value ^= self.toggle_mask;
        }
        value
    }
}

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'T' | 't' => 'A',
            x => x,
        })
        .collect()
}

pub struct BatchSequenceReader {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    dna_db: bool,
    toggle_mask: bool,
    revcom_version: bool,
    str_: Option<String>,
    start_: usize,
    finish_: usize,
    str_pos_: usize,
    queue_: Vec<u64>,
    queue_pos_: usize,
    loaded_ch_: usize,
    last_minimizer_: u64,
    last_ambig_: bool,
    lmer_: u64,
}

impl BatchSequenceReader {
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        dna_db: bool,
        toggle_mask: bool,
        revcom_version: bool,
    ) -> io::Result<Self> {
        Ok(BatchSequenceReader {
            k,
            l,
            spaced_seed_mask,
            dna_db,
            toggle_mask,
            revcom_version,
            str_: None,
            start_: 0,
            finish_: 0,
            str_pos_: 0,
            queue_: Vec::new(),
            queue_pos_: 0,
            loaded_ch_: 0,
            last_minimizer_: !0,
            last_ambig_: false,
            lmer_: 0,
        })
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, len: usize) {
        self.str_ = Some(seq.to_string());
        self.start_ = start;
        self.finish_ = start + len;
        self.str_pos_ = self.start_;
        if ((self.finish_ - self.start_) as isize) < self.l as isize {
            self.str_pos_ = self.finish_;
        }
        self.queue_.clear();
        self.queue_pos_ = 0;
        self.loaded_ch_ = 0;
        self.last_minimizer_ = !0;
        self.last_ambig_ = false;
        self.lmer_ = 0;
    }

    pub fn is_ambiguous(&self) -> bool {
        self.last_ambig_
    }

    pub fn next_minimizer(&mut self) -> Option<u64> {
        // Implementation details here...
        None // Placeholder
    }

    pub fn last_minimizer(&self) -> u64 {
        self.last_minimizer_
    }
}
