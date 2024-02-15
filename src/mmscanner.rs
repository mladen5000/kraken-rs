pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const BITS_PER_CHAR_DNA: usize = 2;
pub const BITS_PER_CHAR_PRO: usize = 4;
pub const CURRENT_REVCOM_VERSION: usize = 1;

pub struct MinimizerScanner {
    lookup_table: [u8; 256],
    k: usize,
    l: usize,
    is_protein: bool,
    mask: u64,
    sequence: Vec<u8>,
    start: usize,
    finish: usize,
}

impl MinimizerScanner {
    pub fn new(k: usize, l: usize, is_protein: bool, mask: u64) -> Self {
        Self {
            lookup_table: [0; 256],
            k,
            l,
            is_protein,
            mask,
            sequence: Vec::new(),
            start: 0,
            finish: 0,
        }
    }

    pub fn set_lookup_table_character(&mut self, c: char, val: u8) {
        self.lookup_table[c as usize] = val;
    }

    pub fn load_sequence(
        &mut self,
        seq: &str,
        start: usize,
        finish: usize,
    ) -> Result<(), &'static str> {
        self.sequence = seq.as_bytes().to_vec();
        self.start = start;
        self.finish = finish;
        Ok(())
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
                let lookup_code = self.lookup_table_[self.str_[self.str_pos_] as usize];
                self.str_pos_ += 1;
                if lookup_code == u8::MAX {
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
                canonical_representation(self.lmer_, self.l_)
            } else {
                self.lmer_
            };
            let mut candidate_lmer = canonical_lmer;
            if let Some(mask) = self.spaced_seed_mask_ {
                candidate_lmer &= mask;
            }
            candidate_lmer ^= self.toggle_mask_;
            if self.k_ == self.l_ as usize {
                self.last_minimizer_ = candidate_lmer ^ self.toggle_mask_;
                return Some(&self.last_minimizer_);
            }
            while let Some(back) = self.queue_.back() {
                if back.candidate > candidate_lmer {
                    self.queue_.pop_back();
                } else {
                    break;
                }
            }
            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos_,
            };
            if self.queue_.is_empty() && self.queue_pos_ >= self.k_ - self.l_ as usize {
                changed_minimizer = true;
            }
            self.queue_.push_back(data);
            if let Some(front) = self.queue_.front() {
                if front.pos < self.queue_pos_ - self.k_ + self.l_ as usize {
                    self.queue_.pop_front();
                    changed_minimizer = true;
                }
            }
            if self.queue_pos_ == self.k_ - self.l_ as usize {
                changed_minimizer = true;
            }
            self.queue_pos_ += 1;
            if self.str_pos_ >= self.k_ {
                break;
            }
        }
        assert!(!self.queue_.is_empty());
        self.last_minimizer_ = self.queue_.front().unwrap().candidate ^ self.toggle_mask_;
        Some(&self.last_minimizer_)
    }

    pub fn reverse_complement(kmer: u64, n: u8) -> u64 {
        let mut rev_comp = 0;
        for i in 0..n {
            let base = (kmer >> (2 * i)) & 3;
            rev_comp |= ((3 - base) << (2 * (n - 1 - i)));
        }
        rev_comp
    }

    pub fn canonical_representation(kmer: u64, n: u8) -> u64 {
        let rev_comp = MinimizerScanner::reverse_complement(kmer, n);
        if kmer < rev_comp {
            kmer
        } else {
            rev_comp
        }
    }

    pub(crate) fn is_ambiguous(&self) -> bool {
        todo!()
    }
}
