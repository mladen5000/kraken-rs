pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const BITS_PER_CHAR_DNA: i32 = 2;
pub const BITS_PER_CHAR_PRO: i32 = 4;
pub const CURRENT_REVCOM_VERSION: i32 = 1;

pub struct MinimizerData {
    candidate: u64,
    pos: isize,
}

pub struct MinimizerScanner {
    k: isize,
    l: isize,
    spaced_seed_mask: u64,
    dna: bool,
    toggle_mask: u64,
    lmer: u64,
    lmer_mask: u64,
    last_minimizer: u64,
    loaded_ch: isize,
    queue: Vec<MinimizerData>,
    queue_pos: isize,
    last_ambig: u64,
    lookup_table: [u8; 256],
    revcom_version: i32,
    str: Option<String>,
    str_pos: usize,
    start: usize,
    finish: usize,
}

impl MinimizerScanner {
    pub fn new(
        k: isize,
        l: isize,
        spaced_seed_mask: u64,
        dna: bool,
        toggle_mask: u64,
        revcom_version: i32,
    ) -> Self {
        let mut scanner = MinimizerScanner {
            str: "",
            k,
            l,
            str_pos: 0,
            start: 0,
            finish: 0,
            spaced_seed_mask,
            dna: dna_sequence,
            toggle_mask,
            loaded_ch: 0,
            last_ambig: 0,
            revcom_version,
            lookup_table: [u8::MAX; 256],
            // Other fields...
        };

        if scanner.l
            > ((std::mem::size_of::<u64>() * 8 - 1)
                / (if scanner.dna {
                    BITS_PER_CHAR_DNA
                } else {
                    BITS_PER_CHAR_PRO
                }))
        {
            panic!(
                "l exceeds size limits for minimizer {} scanner",
                if scanner.dna { "nucleotide" } else { "protein" }
            );
        }

        scanner.lmer_mask = 1;
        scanner.lmer_mask <<= (scanner.l
            * (if scanner.dna {
                BITS_PER_CHAR_DNA
            } else {
                BITS_PER_CHAR_PRO
            }));
        scanner.lmer_mask -= 1;
        scanner.toggle_mask &= scanner.lmer_mask;

        if scanner.finish == usize::MAX {
            scanner.finish = scanner.str.len();
        }

        if (scanner.finish - scanner.start) + 1 < scanner.l as usize {
            scanner.str_pos = scanner.finish;
        }

        if scanner.dna {
            scanner.set_lookup_table_character('A', 0x00);
            scanner.set_lookup_table_character('C', 0x01);
            scanner.set_lookup_table_character('G', 0x02);
            scanner.set_lookup_table_character('T', 0x03);
        } else {
            // Reduced alphabet uses 15-letter alphabet from AD Solis (2015),
            // in Proteins (doi:10.1002/prot.24936)
            // This means that the ambiguous codes B (N/D) and Z (E/Q) aren't
            // usable here because they span multiple groups.  Although J (I/L)
            // could be used because I & L are in the same group, we treat it as
            // we do B, Z, and X for consistency and to ease future changes to
            // the code base.

            // stop codons/rare amino acids
            scanner.set_lookup_table_character('*', 0x00);
            scanner.set_lookup_table_character('U', 0x00);
            scanner.set_lookup_table_character('O', 0x00);
            // alanine
            scanner.set_lookup_table_character('A', 0x01);
            // asparagine, glutamine, serine
            scanner.set_lookup_table_character('N', 0x02);
            scanner.set_lookup_table_character('Q', 0x02);
            scanner.set_lookup_table_character('S', 0x02);
            // cysteine
            scanner.set_lookup_table_character('C', 0x03);
            // aspartic acid, glutamic acid
            scanner.set_lookup_table_character('D', 0x04);
            scanner.set_lookup_table_character('E', 0x04);
            // phenylalanine
            scanner.set_lookup_table_character('F', 0x05);
            // glycine
            scanner.set_lookup_table_character('G', 0x06);
            // histidine
            scanner.set_lookup_table_character('H', 0x07);
            // isoleucine, leucine
            scanner.set_lookup_table_character('I', 0x08);
            scanner.set_lookup_table_character('L', 0x08);
            // lysine
            scanner.set_lookup_table_character('K', 0x09);
            // proline
            scanner.set_lookup_table_character('P', 0x0a);
            // arginine
            scanner.set_lookup_table_character('R', 0x0b);
            // methionine, valine
            scanner.set_lookup_table_character('M', 0x0c);
            scanner.set_lookup_table_character('V', 0x0c);
            // threonine
            scanner.set_lookup_table_character('T', 0x0d);
            // tryptophan
            scanner.set_lookup_table_character('W', 0x0e);
            // tyrosine
            scanner.set_lookup_table_character('Y', 0x0f);
        }

        scanner
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, finish: usize) {
        self.str = seq;
        self.start = start;
        self.finish = finish;
        self.str_pos = start;
        if self.finish > self.str.len() {
            self.finish = self.str.len();
        }
        if (self.finish - self.start) + 1 < self.l as usize {
            self.str_pos = self.finish;
        }
        self.queue.clear();
        self.queue_pos = 0;
        self.loaded_ch = 0;
        self.last_minimizer = !0;
        self.last_ambig = 0;
    }

    pub fn next_minimizer(&mut self) -> Option<&u64> {
        if self.str_pos >= self.finish {
            return None;
        }
        let mut changed_minimizer = false;
        let bits_per_char = if self.dna {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };
        let ambig_code = (1 << bits_per_char) - 1;
        while !changed_minimizer {
            if self.loaded_ch == self.l {
                self.loaded_ch -= 1;
            }
            while self.loaded_ch < self.l && self.str_pos < self.finish {
                self.loaded_ch += 1;
                self.lmer <<= bits_per_char;
                self.last_ambig <<= bits_per_char;
                let lookup_code = lookup_table[(self.str[self.str_pos] as usize)];
                self.str_pos += 1;
                if lookup_code == u8::MAX {
                    self.queue.clear();
                    self.queue_pos = 0;
                    self.lmer = 0;
                    self.loaded_ch = 0;
                    self.last_ambig |= ambig_code;
                } else {
                    self.lmer |= lookup_code;
                }
                self.lmer &= self.lmer_mask;
                self.last_ambig &= self.lmer_mask;
                if (self.str_pos - self.start) >= self.k as usize && self.loaded_ch < self.l {
                    return Some(self.last_minimizer);
                }
            }
            if self.loaded_ch < self.l {
                return None;
            }
            let mut canonical_lmer = if self.dna {
                Self::canonical_representation(self.lmer, self.l)
            } else {
                self.lmer
            };
            if let Some(spaced_seed_mask) = self.spaced_seed_mask {
                canonical_lmer &= spaced_seed_mask;
            }
            let candidate_lmer = canonical_lmer ^ self.toggle_mask;
            if self.k == self.l {
                self.last_minimizer = candidate_lmer ^ self.toggle_mask;
                return Some(self.last_minimizer);
            }
            while let Some(back) = self.queue.back() {
                if back.candidate > candidate_lmer {
                    self.queue.pop_back();
                } else {
                    break;
                }
            }
            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos,
            };
            if self.queue.is_empty() && self.queue_pos >= self.k - self.l {
                changed_minimizer = true;
            }
            self.queue.push_back(data);
            if let Some(front) = self.queue.front() {
                if front.pos < self.queue_pos - self.k + self.l {
                    self.queue.pop_front();
                    changed_minimizer = true;
                }
            }
            if self.queue_pos == self.k - self.l {
                changed_minimizer = true;
            }
            self.queue_pos += 1;
            if self.str_pos >= self.k as usize {
                break;
            }
        }
        assert!(!self.queue.is_empty());
        self.last_minimizer = self.queue.front().unwrap().candidate ^ self.toggle_mask;
        Some(self.last_minimizer)
    }

    pub fn last_minimizer(&self) -> u64 {
        self.last_minimizer
    }

    pub fn k(&self) -> isize {
        self.k
    }

    pub fn l(&self) -> isize {
        self.l
    }

    pub fn is_dna(&self) -> bool {
        self.dna
    }

    pub fn is_ambiguous(&self) -> bool {
        // Implement the method here
        (self.queue_pos < self.k - self.l) || (self.last_ambig != 0)
    }

    pub fn reverse_complement(kmer: u64, n: u8) -> u64 {
        let mut kmer = kmer;
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        if self.revcom_version == 0 {
            return (!kmer) & ((1u64 << (n as u64 * 2)) - 1);
        } else {
            return ((!kmer) >> (std::mem::size_of::<u64>() * 8 - n as usize * 2))
                & ((1u64 << (n as u64 * 2)) - 1);
        }
    }

    pub fn canonical_representation(&self, kmer: u64, n: u8) -> u64 {
        let revcom = self.reverse_complement(kmer, n);
        return if kmer < revcom { kmer } else { revcom };
    }

    fn set_lookup_table_character(&mut self, ch: char, val: u8) {
        // Implement the method here
        self.lookup_table[c as usize] = val;
        self.lookup_table[c.to_ascii_lowercase as usize] = val;
    }
}
