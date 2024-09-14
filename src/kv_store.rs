/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */


pub type HKey = u64;
pub type HValue = u32;

/// A trait representing a key-value store with 64-bit keys and 32-bit values.
pub trait KeyValueStore {
    fn get(&self, key: HKey) -> HValue;
}

/// Implementation of the MurmurHash3 64-bit hash function.
///
/// This function is used to hash 64-bit keys into 64-bit hash codes.
/// The hash function is the 64-bit variant of MurmurHash3.
#[inline]
pub fn murmur_hash3(key: HKey) -> u64 {
    let mut k = key;
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}
