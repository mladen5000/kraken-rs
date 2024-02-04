use std::cmp;
use std::collections::HashSet;
use std::f64;
use std::f64::consts::LN_2;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufWriter, Write};
use std::iter;
use std::num::Wrapping;
use std::ops::{BitAnd, Shl};
use std::path::Path;

// Constants based on the provided C++ code
const P_PRIME: u8 = 25;
const M_PRIME: usize = 1 << P_PRIME;

// Helper functions for bit operations
fn clz(x: u64) -> u8 {
    x.leading_zeros() as u8
}

fn extract_bits<T: Copy + From<u8> + Shl<Output = T> + BitAnd<Output = T>>(
    value: T,
    hi: u8,
    lo: u8,
    shift_left: bool,
) -> T {
    // Implementation similar to the C++ code
    let bitmask: T = (((T::from(1u8) << (hi - lo)) - T::from(1u8)) << lo).into();
    let result = value & bitmask;

    if shift_left {
        result << (std::mem::size_of::<T>() * 8 - hi) as usize
    } else {
        result >> lo as usize
    }
}

// HyperLogLogPlusMinus data structure
pub struct HyperLogLogPlusMinus {
    p: u8,
    m: usize,
    // bit_mixer: fn(u64) -> u64,
    registers: Vec<u8>,
    sparse: bool,
    sparse_list: HashSet<u32>,
    hash_function: fn(u64) -> u64, // Function pointer for the hash function
    pub use_n_observed: bool,
}

impl<T> HyperLogLogPlusMinus<T> {
    // Constructor
    fn new(precision: u8, sparse: bool) -> Self {
        Self {
            p: precision,
            m: 1 << precision,
            sparse,
            registers: vec![0; 1 << precision],
            sparse_list: HashSet::new(),
            hash_function: hash_functions::wang_mixer,
            use_n_observed: true,
            // bit_mixer: Some(default_hash_function),
        }
    }
    pub fn n_observed(&self) -> u64 {
        self.n_observed
    }

    pub fn insert(&mut self, items: Vec<u64>) {
        for i in items {
            self.insert(i);
        }
    }
    pub fn reset(&mut self) {
        self.sparse = true;
        self.sparse_list.clear();
        self.M.clear();
    }
    pub fn merge(&mut self, other: &Self) {
        if self.p != other.p {
            panic!("Cannot merge HyperLogLogPlusMinus instances with different precisions");
        }
        if other.n_observed() == 0 {
            return;
        }
        if self.n_observed() == 0 {
            let n_observed = other.n_observed();
            let sparse = other.sparse;
            let sparselist = other.sparselist;
            self.M = other.M.clone();
        } else {
            let n_observed = self.n_observed + other.n_observed();
            if self.sparse && other.sparse {
                self.sparse_list.extend(other.sparse_list.iter());
            } else if other.sparse {
                self.add_to_registers(other.sparse_list);
            } else {
                if self.sparse {
                    self.sparse = false;
                    self.M = other.M;
                    self.add_to_registers(self.sparse_list);
                    self.sparse_list.clear()
                } else {
                    //merge registers
                    for i in 0..M.size() {
                        if other.M[i] > self.M[i] {
                            self.M[i] = other.M[i];
                        }
                    }
                }
            }
        }
    }
    fn cardinality(&self) -> u64 {
        /* Implementation here */
        ertl_cardinality();
    }

    fn flajolet_cardinality(self, use_sparse_precision: bool) -> () {
        M = vec![0; self.m];
        if use_sparse_precision {
            return linear_counting(self.m_prime, self.m_prime - self.sparse_list.len());
        } else {
            // for testing purposes
            M = vec![0; self.m];
            for val in self.sparse_list {
                let idx = get_index(val, self.p);
                assert!(idx < self.m.len());
                let rank_val = get_encoded_rank(val, self.p_prime, self.p);
                if rank_val > self.M[idx] {
                    self.M[idx] = rank_val;
                }
            }
        }
        let est = calculate_raw_estimate(M);
        if est <= self.m * 2.5 {
            let v = count_zeros(M);
            if v > 0 {
                est = linear_counting(m, v);
            }
        }
        if self.use_n_observed && self.n_observed < est {
            self.n_observed as u64
        } else {
            est.round() as u64
        }
    }

    fn ertl_cardinality(self) -> u64 {
        let (q, m): (usize, usize);
        let C = vec![];

        if self.sparse {
            q = 64 - self.p_prime;
            m = self.m_prime;
            C = sparse_register_histogram(self.sparse_list, self.p_prime, p, q)
        } else {
            q = 64 - self.p;
            m = self.m;
            C = register_histogram(self.m, q);
        }
        let mut est_denominator = m * tau(1.0 - C[q + 1] as f64 / m as f64);
        // skip prints
        // skip prints
        for k in (0..q).rev() {
            est_denominator += C[k] as f64;
            est_denominator *= 0.5;
        }
        // skip prints
        // skip prints
        est_denominator += m * sigma(C[q + 1] as f64 / m as f64);
        // skip prints
        unsafe {
            m_sq_alpha_inf = m / (2.0 * logf64(2)) * m;
        };
        est = m_sq_alpha_inf / est_denominator;

        if self.use_n_observed && self.n_observed < est {
            self.n_observed as u64
        } else {
            est.round() as u64
        }
    }

    fn heule_cardinality(&self, correct_bias: bool) -> u64 {
        if self.p > 18 {
            eprintln!("Huele HLL++ only supports precisions up to 18");
            return self.ertl_cardinality();
        }
        if self.sparse {
            // if sparse then use linear counting
            let lc_estimate = self.linea_counting(
                self.m_prime,
                self.m_prime - self.sparse_list.len(),
                lc_estimate,
            );
            eprintln!(
                "Sparse representation, using linear counting estimate: {}",
                lc_estimate
            );
            return lc_estimate;
        }

        let v = count_zeros(&self.m);
        if v != 0 {
            let lc_esimate = (self.linear_counting(self.m, v) as f64).round() as u64;
            eprintln!("Using linear counting estimate: {}", lc_esimate);
            // check if lc estimate is below the threshold
            // assert(lc_estimate >= 0)
            if lc_esimate <= self.threshold[self.p as usize - 4] {
                return lc_esimate;
            }
            eprintln!(
                "Linear counting estimate above threshold {}",
                self.threshold[self.p as usize - 4]
            );
        }

        // estimate cardinality on registers
        let mut est = calculate_raw_estimate(&self.m);
        eprintln!(" Raw estimate: {}", est);
        // correct for biases if estimate smaller than 5m
        if (correct_bias && est <= self.m * 5.0) {
            let bias = self.get_estimate_bias(est, self.p);
            assert!(est > bias);
            est -= bias;
        }
        if self.use_n_observed && self.n_observed < est {
            self.n_observed as u64
        } else {
            est.round() as u64
        }
    }

    fn cardinality(&self) -> u64 {
        self.ertl_cardinality()
    }

    fn size(&self) -> u64 {
        self.cardinality()
    }

    fn switch_to_normal_representation(&mut self) {
        if !self.sparse {
            return;
        }
        self.sparse = false;
        self.M = vec![0u8; self.m];
        add_to_registers(self.sparse_list);
        self.sparse_list.clear();
    }

    fn add_to_registers(&mut self, sparse_list: SparseListType) {
        if self.sparse {
            eprintln!("Cannot add to registers in sparse representation");
            return;
        }
        if self.sparse_list.size() == 0 {
            return;
        }
        for encoded_hash_value_ptr in (self.sparse_list.begin()..self.sparse_list.end()) {
            let idx: usize = get_index(*encoded_hash_value_ptr, self.p);
            assert(idx < M.size());
            let rank_val = get_encoded_rank(*encoded_hash_value_ptr, P_PRIME, self.p);
            if rank_val > self.M[idx] {
                self.M[idx] = rank_val;
            }
        }
    }
    // Additional methods as needed
}

// Methods for sparse representation, bias correction, and improved estimators
mod ertl_improved_estimator {

    fn sigma(x: f64) -> f64 {
        /* Implementation here */
        assert!(x);
        if x == 1.0 {
            return f64::INFINITY;
        }
        let (prev_sigma_x, sigma_x, y) = (0.0, 0.0, 0.0);

        loop {
            prev_sigma_x = sigma_x;
            x *= x;
            sigma_x += x * y;
            y += y;
            if sigma_x == prev_sigma_x {
                break;
            }
        }
        return sigma_x;
    }

    fn tau(x: f64) -> f64 {
        assert!(x >= 0.0 && x <= 1.0);
        if (x == 0.0 || x == 1.0) {
            return 0.0;
        }
        let (prev_tau_x, y, tau_x) = (0.0, 1.0, 1.0 - x);

        loop {
            prev_tau_x = tau_x;
            x = x.sqrt();
            y /= 2.0;
            tau_x -= std::pow(1 - x, 2) * y;
            if tau_x == prev_tau_x {
                break;
            }
        }
        return tau_x / 3.0;
    }
    // Additional functions as needed
}

// Hash functions and other utilities
pub mod hash_functions {
    pub fn ranhash(u: u64) -> u64 {
        /* Implementation here */
        let mut v: u64 = u * 3935559000370003845 + 2691343689449507681;
        v ^= v >> 21;
        v ^= v << 37;
        v ^= v >> 4;
        v *= 4768777513237032717;
        v ^= v << 20;
        v ^= v >> 41;
        v ^= v << 5;

        return v;
    }
    pub fn murmurhash3_finalizer(mut key: u64) -> u64 {
        //murmurhash returns a hash value of 0 for the key 0 - avoid that.
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccd;
        key ^= key >> 33;
        key *= 0xc4ceb9fe1a85ec53;
        key ^= key >> 33;
        return key;
    }
    pub fn wang_mixer(mut key: u64) -> u64 {
        key = (!key) + (key << 21);
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }
    // Additional hash functions as needed
}

// Reporting functionality

// Define a struct for HyperLogLogPlusMinus
impl HyperLogLogPlusMinus {
    // Constructor with hash function parameter
    fn new(precision: u8, sparse: bool, hash_function: fn(u64) -> u64) -> Self {
        // Initialization, including setting the hash function
    }

    // Methods as defined before, utilizing the hash function
    fn insert(&mut self, item: u64) {
        let hash = (self.hash_function)(item);
        // Insertion logic using the hash value
    }

    // Additional methods as needed
}

// Unit tests for verifying functionality
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_operations() { /* Tests here */
    }
    #[test]
    fn test_hyperloglogplusminus() { /* Tests here */
    }
    // Additional tests as needed
}
