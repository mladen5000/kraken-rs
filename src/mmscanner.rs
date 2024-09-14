/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::VecDeque;

const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
const DEFAULT_SPACED_SEED_MASK: u64 = 0;
const BITS_PER_CHAR_DNA: usize = 2;
const BITS_PER_CHAR_PRO: usize = 4;
const CURRENT_REVCOM_VERSION: i32 = 1;

pub struct MinimizerData {
    candidate: u64,
    pos: isize,
}

pub struct MinimizerScanner<'a> {
    str_: Option<&'a str>, // pointer to sequence
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
    queue_: VecDeque<MinimizerData>,
    queue_pos_: isize,
    last_ambig_: u64,
    lookup_table_: [u8; 256],
    revcom_version_: i32,
}

impl<'a> MinimizerScanner<'a> {
    pub fn new(
        k: isize,
        l: isize,
        spaced_seed_mask: u64,
        dna_sequence: bool,
        toggle_mask: u64,
        revcom_version: i32,
    ) -> Self {
        let mut scanner = MinimizerScanner {
            str_: None,
            k_: k,
            l_: l,
            str_pos_: 0,
            start_: 0,
            finish_: 0,
            spaced_seed_mask_: spaced_seed_mask,
            dna_: dna_sequence,
            toggle_mask_: toggle_mask,
            lmer_: 0,
            lmer_mask_: 0,
            last_minimizer_: u64::MAX,
            loaded_ch_: 0,
            queue_: VecDeque::new(),
            queue_pos_: 0,
            last_ambig_: 0,
            lookup_table_: [u8::MAX; 256],
            revcom_version_: revcom_version,
        };

        let bits_per_char = if scanner.dna_ {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };

        if l > ((std::mem::size_of::<u64>() * 8 - 1) / bits_per_char) as isize {
            panic!(
                "l exceeds size limits for minimizer {} scanner",
                if dna_sequence {
                    "nucleotide"
                } else {
                    "protein"
                }
            );
        }

        scanner.lmer_mask_ = (1 << (l * bits_per_char as isize)) - 1;
        scanner.toggle_mask_ &= scanner.lmer_mask_;

        if scanner.dna_ {
            scanner.set_lookup_table_character('A', 0x00);
            scanner.set_lookup_table_character('C', 0x01);
            scanner.set_lookup_table_character('G', 0x02);
            scanner.set_lookup_table_character('T', 0x03);
        } else {
            scanner.set_lookup_table_character('*', 0x00);
            scanner.set_lookup_table_character('U', 0x00);
            scanner.set_lookup_table_character('O', 0x00);
            scanner.set_lookup_table_character('A', 0x01);
            scanner.set_lookup_table_character('N', 0x02);
            scanner.set_lookup_table_character('Q', 0x02);
            scanner.set_lookup_table_character('S', 0x02);
            scanner.set_lookup_table_character('C', 0x03);
            scanner.set_lookup_table_character('D', 0x04);
            scanner.set_lookup_table_character('E', 0x04);
            scanner.set_lookup_table_character('F', 0x05);
            scanner.set_lookup_table_character('G', 0x06);
            scanner.set_lookup_table_character('H', 0x07);
            scanner.set_lookup_table_character('I', 0x08);
            scanner.set_lookup_table_character('L', 0x08);
            scanner.set_lookup_table_character('K', 0x09);
            scanner.set_lookup_table_character('P', 0x0a);
            scanner.set_lookup_table_character('R', 0x0b);
            scanner.set_lookup_table_character('M', 0x0c);
            scanner.set_lookup_table_character('V', 0x0c);
            scanner.set_lookup_table_character('T', 0x0d);
            scanner.set_lookup_table_character('W', 0x0e);
            scanner.set_lookup_table_character('Y', 0x0f);
        }

        scanner
    }

    pub fn load_sequence(&mut self, seq: &'a str, start: usize, finish: usize) {
        self.str_ = Some(seq);
        self.start_ = start;
        self.finish_ = finish;
        self.str_pos_ = start;
        if self.finish_ > seq.len() {
            self.finish_ = seq.len();
        }
        if (self.finish_ - self.start_) + 1 < self.l_ as usize {
            self.str_pos_ = self.finish_;
        }
        self.queue_.clear();
        self.queue_pos_ = 0;
        self.loaded_ch_ = 0;
        self.last_minimizer_ = u64::MAX;
        self.last_ambig_ = 0;
    }

    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.str_pos_ >= self.finish_ {
            return None;
        }

        let bits_per_char = if self.dna_ {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };

        let ambig_code = (1u32 << bits_per_char) - 1;

        loop {
            let mut changed_minimizer = false;

            if self.loaded_ch_ == self.l_ {
                self.loaded_ch_ -= 1;
            }

            while self.loaded_ch_ < self.l_ && self.str_pos_ < self.finish_ {
                self.loaded_ch_ += 1;
                self.lmer_ <<= bits_per_char;
                self.last_ambig_ <<= bits_per_char;

                let lookup_code =
                    self.lookup_table_[self.str_.unwrap().as_bytes()[self.str_pos_] as usize];
                self.str_pos_ += 1;

                if lookup_code == u8::MAX {
                    self.queue_.clear();
                    self.queue_pos_ = 0;
                    self.lmer_ = 0;
                    self.loaded_ch_ = 0;
                    self.last_ambig_ |= ambig_code as u64;
                } else {
                    self.lmer_ |= lookup_code as u64;
                }

                self.lmer_ &= self.lmer_mask_;
                self.last_ambig_ &= self.lmer_mask_;

                if (self.str_pos_ - self.start_) >= self.k_ as usize && self.loaded_ch_ < self.l_ {
                    return Some(self.last_minimizer_);
                }
            }

            if self.loaded_ch_ < self.l_ {
                return None;
            }

            let mut canonical_lmer = if self.dna_ {
                self.canonical_representation(self.lmer_, self.l_ as u8)
            } else {
                self.lmer_
            };

            if self.spaced_seed_mask_ != 0 {
                canonical_lmer &= self.spaced_seed_mask_;
            }

            let candidate_lmer = canonical_lmer ^ self.toggle_mask_;

            if self.k_ == self.l_ {
                self.last_minimizer_ = candidate_lmer ^ self.toggle_mask_;
                return Some(self.last_minimizer_);
            }

            while !self.queue_.is_empty() && self.queue_.back().unwrap().candidate > candidate_lmer
            {
                self.queue_.pop_back();
            }

            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos_,
            };

            if self.queue_.is_empty() && self.queue_pos_ >= self.k_ - self.l_ {
                changed_minimizer = true;
            }

            self.queue_.push_back(data);

            if self.queue_.front().unwrap().pos < self.queue_pos_ - self.k_ + self.l_ {
                self.queue_.pop_front();
                changed_minimizer = true;
            }

            if self.queue_pos_ == self.k_ - self.l_ {
                changed_minimizer = true;
            }

            self.queue_pos_ += 1;

            if self.str_pos_ >= self.k_ as usize {
                break;
            }

            if changed_minimizer {
                break;
            }
        }

        assert!(!self.queue_.is_empty());
        self.last_minimizer_ = self.queue_.front().unwrap().candidate ^ self.toggle_mask_;
        Some(self.last_minimizer_)
    }

    fn set_lookup_table_character(&mut self, c: char, val: u8) {
        self.lookup_table_[c as usize] = val;
        self.lookup_table_[(c as u8).to_ascii_lowercase() as usize] = val;
    }

    fn reverse_complement(&self, kmer: u64, n: u8) -> u64 {
        let mut kmer = kmer;
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);

        if self.revcom_version_ == 0 {
            (!kmer) & ((1u64 << (n * 2)) - 1)
        } else {
            ((!kmer) >> (64 - n * 2)) & ((1u64 << (n * 2)) - 1)
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
        self.queue_pos_ < self.k_ - self.l_ || self.last_ambig_ != 0
    }
}
