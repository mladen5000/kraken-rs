/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::cmp::min;
use std::usize;

pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const BITS_PER_CHAR_DNA: usize = 2;
pub const BITS_PER_CHAR_PRO: usize = 4;
pub const CURRENT_REVCOM_VERSION: i32 = 1;

pub struct MinimizerData {
    pub candidate: u64,
    pub pos: isize,
}

pub struct MinimizerScanner {
    str_: Option<&'static str>,
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
            last_minimizer_: !0,
            loaded_ch_: 0,
            queue_: Vec::new(),
            queue_pos_: 0,
            last_ambig_: 0,
            lookup_table_: [0xff; 256],
            revcom_version_: revcom_version,
        };

        if l > ((std::mem::size_of::<u64>() * 8 - 1) as isize
            / if dna_sequence {
                BITS_PER_CHAR_DNA as isize
            } else {
                BITS_PER_CHAR_PRO as isize
            }) as isize
        {
            panic!(
                "l exceeds size limits for minimizer {} scanner",
                if dna_sequence {
                    "nucleotide"
                } else {
                    "protein"
                }
            );
        }

        scanner.lmer_mask_ = 1;
        scanner.lmer_mask_ <<= ((l as usize)
            * if dna_sequence {
                BITS_PER_CHAR_DNA
            } else {
                BITS_PER_CHAR_PRO
            }) as u64;
        scanner.lmer_mask_ -= 1;
        scanner.toggle_mask_ &= scanner.lmer_mask_;

        if scanner.finish_ == usize::MAX {
            scanner.finish_ = scanner.str_.map_or(0, |s| s.len());
        }

        if (scanner.finish_ - scanner.start_) + 1 < l as usize {
            scanner.str_pos_ = scanner.finish_;
        }

        if dna_sequence {
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

    fn set_lookup_table_character(&mut self, ch: char, val: u8) {
        self.lookup_table_[ch as usize] = val;
        self.lookup_table_[ch.to_ascii_lowercase() as usize] = val;
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, finish: usize) {
        self.str_ = Some(seq);
        self.start_ = start;
        self.finish_ = finish;
        self.str_pos_ = self.start_;

        if self.finish_ > self.str_.map_or(0, |s| s.len()) {
            self.finish_ = self.str_.map_or(0, |s| s.len());
        }

        if (self.finish_ - self.start_) + 1 < self.l_ as usize {
            self.str_pos_ = self.finish_;
        }

        self.queue_.clear();
        self.queue_pos_ = 0;
        self.loaded_ch_ = 0;
        self.last_minimizer_ = !0;
        self.last_ambig_ = 0;
    }

    pub fn next_minimizer(&mut self) -> Option<&u64> {
        if self.str_pos_ >= self.finish_ {
            return None;
        }

        let mut changed_minimizer = false;
        let bits_per_char = if self.dna_ {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };
        let ambig_code = (1u64 << bits_per_char) - 1;

        while !changed_minimizer {
            if self.loaded_ch_ == self.l_ {
                self.loaded_ch_ -= 1;
            }

            while self.loaded_ch_ < self.l_ && self.str_pos_ < self.finish_ {
                self.loaded_ch_ += 1;
                self.lmer_ <<= bits_per_char;
                self.last_ambig_ <<= bits_per_char;

                let lookup_code = self.lookup_table_
                    [self.str_.as_ref().unwrap().as_bytes()[self.str_pos_] as usize];
                self.str_pos_ += 1;

                if lookup_code == 0xff {
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

                if (self.str_pos_ - self.start_) >= self.k_ as usize && self.loaded_ch_ < self.l_ {
                    return Some(&self.last_minimizer_);
                }
            }

            if self.loaded_ch_ < self.l_ {
                return None;
            }

            let canonical_lmer = if self.dna_ {
                self.canonical_representation(self.lmer_, self.l_ as u8)
            } else {
                self.lmer_
            };

            let candidate_lmer = if self.spaced_seed_mask_ != 0 {
                canonical_lmer & self.spaced_seed_mask_
            } else {
                canonical_lmer ^ self.toggle_mask_
            };

            if self.k_ == self.l_ {
                self.last_minimizer_ = candidate_lmer ^ self.toggle_mask_;
                return Some(&self.last_minimizer_);
            }

            while !self.queue_.is_empty() && self.queue_.last().unwrap().candidate > candidate_lmer
            {
                self.queue_.pop();
            }

            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos_,
            };

            if self.queue_.is_empty() && self.queue_pos_ >= self.k_ - self.l_ {
                changed_minimizer = true;
            }

            self.queue_.push(data);

            if !self.queue_.is_empty()
                && self.queue_.first().unwrap().pos < self.queue_pos_ - self.k_ + self.l_
            {
                self.queue_.remove(0);
                changed_minimizer = true;
            }

            if self.queue_pos_ == self.k_ - self.l_ {
                changed_minimizer = true;
            }

            self.queue_pos_ += 1;

            if self.str_pos_ >= self.k_ as usize {
                break;
            }
        }

        assert!(!self.queue_.is_empty());
        self.last_minimizer_ = self.queue_.first().unwrap().candidate ^ self.toggle_mask_;
        Some(&self.last_minimizer_)
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

    fn reverse_complement(&self, kmer: u64, n: u8) -> u64 {
        let mut kmer = kmer;
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCCu64) >> 2) | ((kmer & 0x3333333333333333u64) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0u64) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0Fu64) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00u64) >> 8) | ((kmer & 0x00FF00FF00FF00FFu64) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000u64) >> 16) | ((kmer & 0x0000FFFF0000FFFFu64) << 16);
        kmer = (kmer >> 32) | (kmer << 32);

        let n_usize = n as usize; // Convert n to usize

        if self.revcom_version_ == 0 {
            (!kmer) & ((1u64 << (n_usize * 2)) - 1) // Use n_usize instead of n
        } else {
            ((!kmer) >> (std::mem::size_of::<u64>() * 8 - n_usize * 2))
                & ((1u64 << (n_usize * 2)) - 1) // Use n_usize instead of n
        }
    }

    fn canonical_representation(&self, kmer: u64, n: u8) -> u64 {
        let revcom = self.reverse_complement(kmer, n);
        min(kmer, revcom)
    }
}
