use anyhow::{bail, Result};
use std::char;
use std::cmp;
use std::string::String;
use std::usize;

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
    str_: Option<String>,
    k_: isize,
    l_: isize,
    str_pos_: usize,
    start_: usize,
    finish_: usize,
    spaced_seed_mask_: u64,
    dna_: bool,
    toggle_mask_: u64,
    lmer_: u64,
    lmer_mask_: u64,
    last_minimizer_: u64,
    loaded_ch_: isize,
    queue_: Vec<MinimizerData>,
    queue_pos_: isize,
    last_ambig_: u64,
    lookup_table_: [u8; 256],
    revcom_version_: i32,
}

impl MinimizerScanner {
    pub fn new(
        k: isize,
        l: isize,
        spaced_seed_mask: u64,
        dna_sequence: bool,
        toggle_mask: u64,
        revcom_version: i32,
    ) -> Result<Self> {
        let bits_per_char = if dna_sequence {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };
        if l > ((std::mem::size_of::<u64>() * 8 - 1) as isize) / (bits_per_char as isize) {
            bail!(
                "l exceeds size limits for minimizer {} scanner",
                if dna_sequence {
                    "nucleotide"
                } else {
                    "protein"
                }
            );
        }

        let mut scanner = MinimizerScanner {
            str_: None,
            k_: k,
            l_: l,
            str_pos_: 0,
            start_: 0,
            finish_: usize::MAX,
            spaced_seed_mask_: spaced_seed_mask,
            dna_: dna_sequence,
            toggle_mask_: toggle_mask,
            lmer_: 0,
            lmer_mask_: 0,
            last_minimizer_: !0,
            loaded_ch_: 0,
            queue_: Vec::new(),
            queue_pos_: 0,
            last_ambig_: 0,
            lookup_table_: [u8::MAX; 256],
            revcom_version_: revcom_version,
        };

        scanner.lmer_mask_ = 1;
        scanner.lmer_mask_ <<= (l * (bits_per_char as isize)) as u64;
        scanner.lmer_mask_ -= 1;
        scanner.toggle_mask_ &= scanner.lmer_mask_;

        // Set lookup table defaults
        for i in 0..=u8::MAX {
            scanner.lookup_table_[i as usize] = u8::MAX;
        }

        if dna_sequence {
            scanner.set_lookup_table_character('A', 0x00);
            scanner.set_lookup_table_character('C', 0x01);
            scanner.set_lookup_table_character('G', 0x02);
            scanner.set_lookup_table_character('T', 0x03);
        } else {
            // protein reduced alphabet
            // according to the mapping given in mmscanner.cc
            scanner.set_lookup_table_character('*', 0x00);
            scanner.set_lookup_table_character('U', 0x00);
            scanner.set_lookup_table_character('O', 0x00);
            scanner.set_lookup_table_character('A', 0x01);
            // N,Q,S = 0x02
            scanner.set_lookup_table_character('N', 0x02);
            scanner.set_lookup_table_character('Q', 0x02);
            scanner.set_lookup_table_character('S', 0x02);
            // C = 0x03
            scanner.set_lookup_table_character('C', 0x03);
            // D,E = 0x04
            scanner.set_lookup_table_character('D', 0x04);
            scanner.set_lookup_table_character('E', 0x04);
            // F = 0x05
            scanner.set_lookup_table_character('F', 0x05);
            // G = 0x06
            scanner.set_lookup_table_character('G', 0x06);
            // H = 0x07
            scanner.set_lookup_table_character('H', 0x07);
            // I,L = 0x08
            scanner.set_lookup_table_character('I', 0x08);
            scanner.set_lookup_table_character('L', 0x08);
            // K = 0x09
            scanner.set_lookup_table_character('K', 0x09);
            // P = 0x0a
            scanner.set_lookup_table_character('P', 0x0a);
            // R = 0x0b
            scanner.set_lookup_table_character('R', 0x0b);
            // M,V = 0x0c
            scanner.set_lookup_table_character('M', 0x0c);
            scanner.set_lookup_table_character('V', 0x0c);
            // T = 0x0d
            scanner.set_lookup_table_character('T', 0x0d);
            // W = 0x0e
            scanner.set_lookup_table_character('W', 0x0e);
            // Y = 0x0f
            scanner.set_lookup_table_character('Y', 0x0f);
        }

        Ok(scanner)
    }

    fn set_lookup_table_character(&mut self, ch: char, val: u8) {
        let upper = ch as u8;
        let lower = (ch.to_ascii_lowercase()) as u8;
        self.lookup_table_[upper as usize] = val;
        self.lookup_table_[lower as usize] = val;
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, finish: usize) {
        self.str_ = Some(seq.to_string());
        self.start_ = start;
        self.finish_ = if finish == usize::MAX {
            seq.len()
        } else {
            cmp::min(finish, seq.len())
        };
        self.str_pos_ = self.start_;
        if ((self.finish_ - self.start_) as isize) + 1 < self.l_ {
            // Invalidate scanner if interval < 1 l-mer
            self.str_pos_ = self.finish_;
        }
        self.queue_.clear();
        self.queue_pos_ = 0;
        self.loaded_ch_ = 0;
        self.last_minimizer_ = !0;
        self.last_ambig_ = 0;
        self.lmer_ = 0;
    }

    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.str_.is_none() {
            return None;
        }
        let seq = self.str_.as_ref().unwrap();
        if self.str_pos_ >= self.finish_ {
            // exhausted string
            return None;
        }
        let bits_per_char = if self.dna_ {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };
        let ambig_code = (1u64 << bits_per_char) - 1;
        let mut changed_minimizer = false;

        while !changed_minimizer {
            // Incorporate next char and fill l-mer
            if self.loaded_ch_ == self.l_ {
                self.loaded_ch_ -= 1;
            }
            while self.loaded_ch_ < self.l_ && self.str_pos_ < self.finish_ {
                self.loaded_ch_ += 1;
                self.lmer_ <<= bits_per_char;
                self.last_ambig_ <<= bits_per_char;

                let c = seq.as_bytes()[self.str_pos_];
                self.str_pos_ += 1;
                let lookup_code = self.lookup_table_[c as usize];
                if lookup_code == u8::MAX {
                    // ambiguous char
                    self.queue_.clear();
                    self.queue_pos_ = 0;
                    self.lmer_ = 0;
                    self.loaded_ch_ = 0;
                    self.last_ambig_ |= ambig_code;
                } else {
                    self.lmer_ |= lookup_code as u64;
                }
                self.lmer_ &= self.lmer_mask_;
                self.last_ambig_ &= self.lmer_mask_;

                if ((self.str_pos_ - self.start_) as isize) >= self.k_ && self.loaded_ch_ < self.l_
                {
                    return Some(self.last_minimizer_);
                }
            }

            if self.loaded_ch_ < self.l_ {
                // can't fill l-mer, end
                return None;
            }

            let canonical_lmer = if self.dna_ {
                self.canonical_representation(self.lmer_, self.l_ as u8)
            } else {
                self.lmer_
            };
            let mut candidate_lmer = canonical_lmer;
            if self.spaced_seed_mask_ != 0 {
                candidate_lmer &= self.spaced_seed_mask_;
            }
            candidate_lmer ^= self.toggle_mask_;

            if self.k_ == self.l_ {
                // short-circuit queue work
                self.last_minimizer_ = candidate_lmer ^ self.toggle_mask_;
                return Some(self.last_minimizer_);
            }

            // Sliding window minimum
            while !self.queue_.is_empty() && self.queue_.last().unwrap().candidate > candidate_lmer
            {
                self.queue_.pop();
            }
            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos_,
            };
            let was_empty = self.queue_.is_empty();
            self.queue_.push(data);

            // expire l-mer not in current window
            if !self.queue_.is_empty() && self.queue_[0].pos < self.queue_pos_ - self.k_ + self.l_ {
                self.queue_.remove(0);
                changed_minimizer = true;
            }

            if was_empty && self.queue_pos_ >= self.k_ - self.l_ {
                changed_minimizer = true;
            }

            if self.queue_pos_ == self.k_ - self.l_ {
                changed_minimizer = true;
            }

            self.queue_pos_ += 1;

            if (self.str_pos_ as isize) >= self.k_ {
                // we can return only after first k-mer loaded
                break;
            }
        }

        if self.queue_.is_empty() {
            return None;
        }

        self.last_minimizer_ = self.queue_[0].candidate ^ self.toggle_mask_;
        Some(self.last_minimizer_)
    }

    pub fn last_minimizer(&self) -> u64 {
        self.last_minimizer_
    }

    pub fn k(&self) -> isize {
        self.k_
    }

    pub fn l(&self) -> isize {
        self.l_
    }

    pub fn is_dna(&self) -> bool {
        self.dna_
    }

    pub fn is_ambiguous(&self) -> bool {
        (self.queue_pos_ < self.k_ - self.l_) || (self.last_ambig_ != 0)
    }

    fn reverse_complement(&self, kmer: u64, n: u8) -> u64 {
        let mut kmer = kmer;
        // reverse bits 2 by 2
        // apply transformations as in original code
        // This code is the same bit-twiddling from original:
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        if self.revcom_version_ == 0 {
            ((!kmer) & ((1u64 << (n as u64 * 2)) - 1))
        } else {
            (((!kmer) >> (std::mem::size_of::<u64>() * 8 - n as usize * 2))
                & ((1u64 << (n as u64 * 2)) - 1))
        }
    }

    fn canonical_representation(&self, kmer: u64, n: u8) -> u64 {
        let revcom = self.reverse_complement(kmer, n);
        if kmer < revcom {
            kmer
        } else {
            revcom
        }
    }
}
