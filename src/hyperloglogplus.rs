// hyperloglogplus.rs

use std::collections::HashSet;



/// HyperLogLog++ class for counting the number of unique 64-bit values in a stream.
pub struct HyperLogLogPlusMinus {
    p: u8,                         // precision, set in constructor
    m: usize,                      // number of registers (2^p)
    registers: Vec<u8>,            // registers, size m
    n_observed: u64,               // number of observed items
    sparse: bool,                  // sparse representation of the data?
    pub sparse_list: HashSet<u32>, // sparse representation using a HashSet
    bit_mixer: fn(u64) -> u64,     // hash function for mixing bits
    use_n_observed: bool,          // return min(estimate, n_observed) instead of estimate
}

impl HyperLogLogPlusMinus {
    const P_PRIME: u8 = 25; // precision when using a sparse representation
    const M_PRIME: u32 = 1 << 25; // 2^P_PRIME

    /// Constructs a new HyperLogLog++ instance with the specified precision.
    ///
    /// # Arguments
    ///
    /// * `precision` - The precision of the HyperLogLog++ instance (default: 12).
    /// * `sparse` - Whether to use sparse representation (default: true).
    /// * `bit_mixer` - The hash function for mixing bits (default: murmurhash3_finalizer).
    ///
    /// # Example
    ///
    /// ```
    /// let hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// ```
    pub fn new(precision: u8, sparse: bool, bit_mixer: fn(u64) -> u64) -> Self {
        assert!(
            precision >= 4 && precision <= 18,
            "precision must be between 4 and 18"
        );
        let m = 1 << precision;
        HyperLogLogPlusMinus {
            p: precision,
            m,
            registers: vec![0; m],
            n_observed: 0,
            sparse,
            sparse_list: HashSet::new(),
            bit_mixer,
            use_n_observed: true,
        }
    }

    /// Resets the HyperLogLog++ instance to its initial state.
    pub fn reset(&mut self) {
        self.sparse = true;
        self.sparse_list.clear();
        self.registers.clear();
    }

    /// Inserts a single item into the HyperLogLog++ instance.
    ///
    /// # Arguments
    ///
    /// * `item` - The item to be inserted.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll.insert(42);
    /// ```
    pub fn insert(&mut self, item: u64) {
        self.n_observed += 1;
        let hash_value = (self.bit_mixer)(item);

        if self.sparse && (self.sparse_list.len() + 1) > (self.m / 4) {
            self.switch_to_normal_representation();
        }

        if self.sparse {
            let encoded_hash_value = self.encode_hash_in_32bit(hash_value);
            self.add_hash_to_sparse_list(encoded_hash_value);
        } else {
            let idx = self.get_index(hash_value);
            let rank = self.get_rank(hash_value);

            if rank > self.registers[idx] {
                self.registers[idx] = rank;
            }
        }
    }

    /// Inserts a vector of items into the HyperLogLog++ instance.
    ///
    /// # Arguments
    ///
    /// * `items` - The vector of items to be inserted.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// let items = vec![1, 2, 3, 4, 5];
    /// hll.insert_vec(&items);
    /// ```
    pub fn insert_vec(&mut self, items: &[u64]) {
        for item in items {
            self.insert(*item);
        }
    }

    /// Merges another HyperLogLog++ instance into the current instance.
    ///
    /// # Arguments
    ///
    /// * `other` - The other HyperLogLog++ instance to be merged.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll1 = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// let mut hll2 = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll1.insert(1);
    /// hll2.insert(2);
    /// hll1.merge(&hll2);
    /// ```
    pub fn merge(&mut self, other: &HyperLogLogPlusMinus) {
        assert_eq!(self.p, other.p, "precisions must be equal");
        if other.n_observed == 0 {
            return;
        }

        if self.n_observed == 0 {
            self.n_observed = other.n_observed;
            self.sparse = other.sparse;
            self.sparse_list = other.sparse_list.clone();
            self.registers = other.registers.clone();
        } else {
            self.n_observed += other.n_observed;
            if self.sparse && other.sparse {
                self.sparse_list.extend(other.sparse_list.iter().cloned());
            } else if other.sparse {
                self.add_to_registers(&other.sparse_list);
            } else {
                if self.sparse {
                    self.sparse = false;
                    self.registers = other.registers.clone();

                    // Create a temporary variable to hold the sparse list
                    let sparse_list = std::mem::take(&mut self.sparse_list);

                    self.add_to_registers(&sparse_list);
                } else {
                    for i in 0..other.registers.len() {
                        if other.registers[i] > self.registers[i] {
                            self.registers[i] = other.registers[i];
                        }
                    }
                }
            }
        }
    }

    /// Returns the cardinality estimate using the Ertl estimator.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll.insert(1);
    /// hll.insert(2);
    /// hll.insert(3);
    /// let cardinality = hll.cardinality();
    /// ```
    pub fn cardinality(&self) -> u64 {
        self.ertl_cardinality()
    }

    /// Returns the cardinality estimate using the Heule estimator.
    ///
    /// # Arguments
    ///
    /// * `correct_bias` - Whether to apply bias correction (default: true).
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll.insert(1);
    /// hll.insert(2);
    /// hll.insert(3);
    /// let cardinality = hll.heule_cardinality(true);
    /// ```
    pub fn heule_cardinality(&self, _correct_bias: bool) -> u64 {
        // Implementation of the Heule estimator
        // ...
        unimplemented!()
    }

    /// Returns the cardinality estimate using the Ertl estimator.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll.insert(1);
    /// hll.insert(2);
    /// hll.insert(3);
    /// let cardinality = hll.ertl_cardinality();
    /// ```
    pub fn ertl_cardinality(&self) -> u64 {
        // Implementation of the Ertl estimator
        // ...
        unimplemented!()
    }

    /// Returns the cardinality estimate using the Flajolet estimator.
    ///
    /// # Arguments
    ///
    /// * `use_sparse_precision` - Whether to use sparse precision (default: true).
    ///
    /// # Example
    ///
    /// ```
    /// let mut hll = HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer);
    /// hll.insert(1);
    /// hll.insert(2);
    /// hll.insert(3);
    /// let cardinality = hll.flajolet_cardinality(true);
    /// ```
    pub fn flajolet_cardinality(&self, _use_sparse_precision: bool) -> u64 {
        // Implementation of the Flajolet estimator
        // ...
        unimplemented!()
    }

    /// Returns the number of observed items.
    pub fn n_observed(&self) -> u64 {
        self.n_observed
    }

    // Private methods
    // ...
}

// Hash functions and other utility functions
// ...

// Private methods

impl HyperLogLogPlusMinus {
    fn switch_to_normal_representation(&mut self) {
        if !self.sparse {
            return;
        }
        self.sparse = false;
        self.registers = vec![0; self.m];
        let sparse_list = std::mem::take(&mut self.sparse_list);
        self.add_to_registers(&sparse_list);
        // self.sparse_list.clear();
    }

    fn add_to_registers(&mut self, sparse_list: &HashSet<u32>) {
        if self.sparse {
            panic!("Cannot add to registers of a sparse HLL");
        }
        for &encoded_hash_value in sparse_list {
            let idx = self.get_index_from_encoded(encoded_hash_value);
            let rank = self.get_encoded_rank(encoded_hash_value);
            if rank > self.registers[idx] {
                self.registers[idx] = rank;
            }
        }
    }

    fn encode_hash_in_32bit(&self, hash_value: u64) -> u32 {
        let idx = self.get_index(hash_value) as u32;
        let idx_shifted = idx << (32 - Self::P_PRIME);

        if idx_shifted << self.p == 0 {
            let additional_rank = self.get_rank(hash_value << Self::P_PRIME);
            idx_shifted | ((additional_rank as u32) << 1) | 1
        } else {
            idx_shifted
        }
    }

    fn add_hash_to_sparse_list(&mut self, encoded_hash_value: u32) {
        self.sparse_list.insert(encoded_hash_value);
    }

    fn get_index(&self, hash_value: u64) -> usize {
        (hash_value >> (64 - self.p)) as usize
    }

    fn get_index_from_encoded(&self, encoded_hash_value: u32) -> usize {
        (encoded_hash_value >> (32 - self.p)) as usize
    }

    fn get_rank(&self, hash_value: u64) -> u8 {
        let rank_bits = hash_value << self.p;
        self.count_leading_zeros(rank_bits, 64 - self.p) + 1
    }

    fn get_encoded_rank(&self, encoded_hash_value: u32) -> u8 {
        if encoded_hash_value & 1 == 1 {
            let additional_rank = Self::P_PRIME - self.p;
            additional_rank + ((encoded_hash_value >> 1) & 0x3F) as u8
        } else {
            self.get_rank((encoded_hash_value as u64) << 32)
        }
    }

    fn count_leading_zeros(&self, value: u64, max: u8) -> u8 {
        if value == 0 {
            max
        } else {
            value.leading_zeros() as u8
        }
    }
}

// Hash functions and other utility functions

pub fn murmurhash3_finalizer(key: u64) -> u64 {
    let mut key = key.wrapping_add(1);
    key ^= key >> 33;
    key = key.wrapping_mul(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key = key.wrapping_mul(0xc4ceb9fe1a85ec53);
    key ^= key >> 33;
    key
}

fn wang_mixer(key: u64) -> u64 {
    let mut key = !key + (key << 21);
    key ^= key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key ^= key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key ^= key >> 28;
    key + (key << 31)
}
