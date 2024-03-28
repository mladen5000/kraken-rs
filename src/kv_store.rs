// File: kraken2_kv_store.rs/*
/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

pub type HKey = u64;
pub type HValue = u32;

pub trait KeyValueStore {
    fn get(&self, key: HKey) -> HValue;
}

pub fn murmur_hash3(key: HKey) -> u64 {
    let mut k = key;
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}
