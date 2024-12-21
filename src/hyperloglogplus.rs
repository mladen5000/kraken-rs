/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file was originally developed for the KrakenUniq taxonomic classification system.
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

/*
 * Implementation of 64-bit HyperLogLog algorithm by Flajolet et al.,
 * with sparse mode for increased precision with low cardinalities (Stefan Heule et al.), and
 * an improved estimator that does not rely on empirical bias correction data (Otmar Ertl)
 */

use std::collections::{BTreeSet, HashMap};
use std::hash::BuildHasherDefault;
use std::io;
use std::marker::PhantomData;

use std::cmp::Ordering;
use std::f64::consts::LN_2;
use std::mem;
use std::u32;

use std::fmt::Debug;

/// Type alias for the hash function
pub type HashFunction = fn(u64) -> u64;

/// SparseListType can be a BTreeSet for sorted elements
type SparseListType = BTreeSet<u32>;

/// HyperLogLogPlusMinus class for counting the number of unique 64-bit values in a stream
/// Note that only HASH = u64 is implemented.
pub struct HyperLogLogPlusMinus<HASH = u64> {
    p: u8,           // Precision
    m: usize,        // Number of registers
    M: Vec<u8>,      // Registers
    n_observed: u64, // Number of observed elements

    sparse: bool,                // Sparse representation flag
    sparse_list: SparseListType, // Sparse list for low cardinalities

    // Sparse versions of p and m
    p_prime: u8,  // Precision when using sparse representation
    m_prime: u32, // Number of registers in sparse mode

    bit_mixer: HashFunction, // Hash function
    use_n_observed: bool,    // Return min(estimate, n_observed)

    _phantom: PhantomData<HASH>, // To keep the type parameter
}

impl HyperLogLogPlusMinus<u64> {
    /// Constructs a new HyperLogLogPlusMinus with the given precision and optional hash function.
    pub fn new(precision: u8, sparse: bool, bit_mixer: Option<HashFunction>) -> io::Result<Self> {
        if precision > 18 || precision < 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "precision must be between 4 and 18",
            ));
        }

        let m = 1 << precision;
        let p_prime = 25;
        let m_prime = 1 << p_prime;

        Ok(Self {
            p: precision,
            m,
            M: if sparse { Vec::new() } else { vec![0; m] },
            n_observed: 0,
            sparse,
            sparse_list: SparseListType::new(),
            p_prime,
            m_prime,
            bit_mixer: bit_mixer.unwrap_or(murmurhash3_finalizer),
            use_n_observed: true,
            _phantom: PhantomData,
        })
    }

    /// Resets the HyperLogLogPlusMinus to its initial state.
    pub fn reset(&mut self) {
        self.sparse = true;
        self.n_observed = 0;
        self.sparse_list.clear();
        self.M.clear();
    }

    /// Inserts an item into the HyperLogLogPlusMinus.
    pub fn insert(&mut self, item: u64) {
        self.n_observed += 1;
        let hash_value = (self.bit_mixer)(item);

        if self.sparse && self.sparse_list.len() + 1 > self.m / 4 {
            self.switch_to_normal_representation();
        }

        if self.sparse {
            let encoded_hash_value = encode_hash_in_32bit(hash_value, self.p_prime, self.p);
            self.add_hash_to_sparse_list(encoded_hash_value);
        } else {
            let idx = get_index(hash_value, self.p);
            let rank = get_rank(hash_value, self.p);
            if rank > self.M[idx] {
                self.M[idx] = rank;
            }
        }
    }

    /// Inserts a vector of items into the HyperLogLogPlusMinus.
    pub fn insert_multiple(&mut self, items: &[u64]) {
        for &item in items {
            self.insert(item);
        }
    }

    /// Merges another HyperLogLogPlusMinus into this one.
    pub fn merge(&mut self, other: &HyperLogLogPlusMinus<u64>) -> io::Result<()> {
        if self.p != other.p {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "precisions must be equal",
            ));
        }

        if other.n_observed == 0 {
            return Ok(());
        }

        if self.n_observed == 0 {
            self.n_observed = other.n_observed;
            self.sparse = other.sparse;
            self.sparse_list = other.sparse_list.clone();
            self.M = other.M.clone();
        } else {
            self.n_observed += other.n_observed;
            if self.sparse && other.sparse {
                self.sparse_list.extend(&other.sparse_list);
                if self.sparse_list.len() + 1 > self.m / 4 {
                    self.switch_to_normal_representation();
                }
            } else if other.sparse {
                self.add_to_registers(&other.sparse_list);
            } else {
                if self.sparse {
                    self.switch_to_normal_representation();
                }
                for (i, &val) in other.M.iter().enumerate() {
                    if val > self.M[i] {
                        self.M[i] = val;
                    }
                }
            }
        }
        Ok(())
    }

    /// Returns the cardinality estimate using Ertl's method.
    pub fn cardinality(&self) -> u64 {
        self.ertl_cardinality()
    }

    /// Returns the number of observed items.
    pub fn n_observed(&self) -> u64 {
        self.n_observed
    }

    /// Returns the estimated size.
    pub fn size(&self) -> u64 {
        self.cardinality()
    }

    // Private helper methods

    fn switch_to_normal_representation(&mut self) {
        if !self.sparse {
            return;
        }
        self.sparse = false;
        self.M = vec![0; self.m];
        let sparse_list = std::mem::take(&mut self.sparse_list);
        self.add_to_registers(&sparse_list);
    }

    fn add_to_registers(&mut self, sparse_list: &SparseListType) {
        for &encoded_hash_value in sparse_list {
            let idx = get_index_u32(encoded_hash_value, self.p);
            let rank_val = get_encoded_rank(encoded_hash_value, self.p_prime, self.p);
            if rank_val > self.M[idx] {
                self.M[idx] = rank_val;
            }
        }
    }

    fn add_hash_to_sparse_list(&mut self, val: u32) {
        self.sparse_list.insert(val);
    }

    fn ertl_cardinality(&self) -> u64 {
        let (q, m, C) = if self.sparse {
            let q = 64 - self.p_prime;
            let m = self.m_prime as usize;
            let C = sparse_register_histogram(&self.sparse_list, self.p_prime, self.p, q);
            (q, m, C)
        } else {
            let q = 64 - self.p;
            let m = self.m;
            let C = register_histogram(&self.M, q);
            (q, m, C)
        };

        let mut est_denominator = m as f64 * tau(1.0 - (C[q as usize + 1] as f64) / m as f64);
        for k in (1..=q).rev() {
            est_denominator += C[k as usize] as f64;
            est_denominator *= 0.5;
        }
        est_denominator += m as f64 * sigma(C[0] as f64 / m as f64);
        let m_sq_alpha_inf = (m as f64 / (2.0 * LN_2)) * m as f64;
        let est = m_sq_alpha_inf / est_denominator;

        if self.use_n_observed && self.n_observed < est as u64 {
            self.n_observed
        } else {
            est.round() as u64
        }
    }
}

// Helper functions

/// Counts the number of leading zeros in a u64 value.
fn clz(x: u64) -> u8 {
    x.leading_zeros() as u8
}

/// Gets the index from the hash value.
fn get_index(hash_value: u64, p: u8) -> usize {
    (hash_value >> (64 - p)) as usize
}

/// Gets the index from the encoded hash value (u32).
fn get_index_u32(hash_value: u32, p: u8) -> usize {
    (hash_value >> (32 - p)) as usize
}

/// Gets the rank from the hash value.
fn get_rank(hash_value: u64, p: u8) -> u8 {
    let shifted = hash_value << p;
    let rank = clz(shifted) + 1;
    rank
}

/// Gets the encoded rank from the encoded hash value.
fn get_encoded_rank(encoded_hash_value: u32, p_prime: u8, p: u8) -> u8 {
    if (encoded_hash_value & 1) == 1 {
        let additional_rank = p_prime - p;
        let rank = ((encoded_hash_value >> 1) & 0x3F) as u8; // 6 bits for rank
        additional_rank + rank
    } else {
        get_rank_u32(encoded_hash_value, p)
    }
}

/// Gets the rank from a u32 hash value.
fn get_rank_u32(hash_value: u32, p: u8) -> u8 {
    let shifted = hash_value << p;
    let rank = clz(shifted as u64) as u8 - 31 + p;
    rank
}

/// Encodes the hash value into 32 bits for the sparse representation.
fn encode_hash_in_32bit(hash_value: u64, p_prime: u8, p: u8) -> u32 {
    let idx = (hash_value >> (64 - p_prime)) as u32;

    if (idx << p) == 0 {
        let additional_rank = get_rank(hash_value, p_prime);
        idx | ((additional_rank as u32) << 1) | 1
    } else {
        idx
    }
}

/// Computes sigma(x) as per Ertl's method.
fn sigma(mut x: f64) -> f64 {
    if x == 1.0 {
        return f64::INFINITY;
    }

    let mut prev_sigma_x;
    let mut sigma_x = x;
    let mut y = 1.0;

    loop {
        prev_sigma_x = sigma_x;
        x *= x;
        sigma_x += x * y;
        y += y;
        if sigma_x == prev_sigma_x {
            break;
        }
    }
    sigma_x
}

/// Computes tau(x) as per Ertl's method.
fn tau(mut x: f64) -> f64 {
    if x == 0.0 || x == 1.0 {
        return 0.0;
    }

    let mut prev_tau_x;
    let mut tau_x = 1.0 - x;
    let mut y = 1.0;

    loop {
        prev_tau_x = tau_x;
        x = x.sqrt();
        y /= 2.0;
        tau_x -= (1.0 - x).powi(2) * y;
        if tau_x == prev_tau_x {
            break;
        }
    }
    tau_x / 3.0
}

/// Generates the register histogram for Ertl's method.
fn register_histogram(M: &[u8], q: u8) -> Vec<u32> {
    let mut C = vec![0u32; q as usize + 2];
    for &val in M {
        if (val as usize) >= C.len() {
            continue;
        }
        C[val as usize] += 1;
    }
    C
}

/// Generates the sparse register histogram for Ertl's method.
fn sparse_register_histogram(sparse_list: &SparseListType, p_prime: u8, p: u8, q: u8) -> Vec<u32> {
    let mut C = vec![0u32; q as usize + 2];
    let mut m = 1u32 << p_prime;
    for &encoded_hash_value in sparse_list {
        let rank_val = get_encoded_rank(encoded_hash_value, p_prime, p);
        if (rank_val as usize) >= C.len() {
            continue;
        }
        C[rank_val as usize] += 1;
        m -= 1;
    }
    C[0] = m;
    C
}

// Hash functions

/// MurmurHash3 finalizer for 64-bit inputs.
pub fn murmurhash3_finalizer(mut key: u64) -> u64 {
    key = key.wrapping_add(1);
    key ^= key >> 33;
    key = key.wrapping_mul(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key = key.wrapping_mul(0xc4ceb9fe1a85ec53);
    key ^= key >> 33;
    key
}

/// Wang's 64-bit hash function.
pub fn wang_mixer(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8);
    key ^= key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4);
    key ^= key >> 28;
    key = key.wrapping_add(key << 31);
    key
}
