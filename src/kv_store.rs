pub struct KeyValueStore {}

pub type HKey = u64;
pub type HValue = u32;

impl KeyValueStore {
    pub fn get(&self, key: HKey) -> Option<HValue> {
        Some(HValue::from(0))
    }
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
