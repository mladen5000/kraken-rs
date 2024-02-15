// File: kraken2_kv_store.rs

use crate::kraken2_data;
use std::hash::{Hash, Hasher};
// Type aliases
type HKey = u64;
type HValue = u32;

// The Key-Value Store Trait
trait KeyValueStore {
    fn get(&self, key: HKey) -> HValue;
}

// Rust-idiomatic hash function, consider using DefaultHasher if suitable
fn murmur_hash3(key: HKey) -> u64 {
    let mut hasher = MurmurHasher::default();
    key.hash(&mut hasher);
    hasher.finish()
}

// A custom MurmurHash3 implementation as in the original C++
struct MurmurHasher {
    state: u64,
}

impl Default for MurmurHasher {
    fn default() -> Self {
        MurmurHasher { state: 0 }
    }
}

impl Hasher for MurmurHasher {
    fn finish(&self) -> u64 {
        self.state
    }

    fn write(&mut self, bytes: &[u8]) {
        let k = u64::from_le_bytes(bytes.try_into().unwrap());

        self.state ^= k >> 33;
        self.state *= 0xff51afd7ed558ccd;
        self.state ^= self.state >> 33;
        self.state *= 0xc4ceb9fe1a85ec53;
        self.state ^= self.state >> 33;
    }
}
