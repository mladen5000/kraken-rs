struct KeyValueStore {}

pub type HKey = u64;
pub type HValue = u32;

pub trait KeyValueStore {
    fn get(&self, key: HKey) -> Option<HValue>;
}

pub fn murmurhash3(key: u64) -> u64 {
    let mut k: u64 = key as u64;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    k
}
