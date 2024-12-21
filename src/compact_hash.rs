use anyhow::{bail, Result};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::sync::Mutex;

// Type aliases, inferred from code usage:
type hkey_t = u64;
type hvalue_t = u64;

// Constants from the code (not given explicitly, we assume from code)
const LOCK_ZONES: usize = 64; // number of lock zones, arbitrary guess
                              // If code used #ifdef LINEAR_PROBING, we can define a feature to toggle.
#[cfg(feature = "linear_probing")]
const LINEAR_PROBING: bool = true;
#[cfg(not(feature = "linear_probing"))]
const LINEAR_PROBING: bool = false;

// A compact hash cell stores key and value in a single u64.
// The value occupies the lower `value_bits` bits, the key occupies the upper `key_bits` bits.
#[derive(Copy, Clone)]
struct CompactHashCell {
    data: u64,
}

impl CompactHashCell {
    fn new() -> Self {
        CompactHashCell { data: 0 }
    }

    fn value(&self, value_bits: usize) -> hvalue_t {
        let mask = (1u64 << value_bits) - 1;
        self.data & mask
    }

    fn hashed_key(&self, value_bits: usize) -> u64 {
        self.data >> value_bits
    }

    fn populate(
        &mut self,
        compacted_key: u64,
        new_value: hvalue_t,
        key_bits: usize,
        value_bits: usize,
    ) {
        // data = (compacted_key << value_bits) | new_value
        self.data = (compacted_key << value_bits) | (new_value & ((1 << value_bits) - 1));
    }
}

pub struct CompactHashTable {
    capacity_: usize,
    size_: usize,
    key_bits_: usize,
    value_bits_: usize,
    file_backed_: bool,
    locks_initialized_: bool,
    table_: *mut CompactHashCell,
    zone_locks_: Vec<Mutex<()>>,
    // If memory mapping is desired, store file and mapped memory here.
    // For now, we emulate with normal load/store.
    backing_data_: Vec<CompactHashCell>,
}

// taxon_counts_t was a map of value -> count
type taxon_counts_t = HashMap<hvalue_t, u64>;

// Mock MurmurHash3 function. The original code uses MurmurHash3(key) returning u64.
// Implementing a real MurmurHash3 is possible but omitted here.
fn murmur_hash3(key: u64) -> u64 {
    // A placeholder hash. In production, implement properly.
    // For stable hashing, consider a standard hash function.
    // But to mimic original code, just do something simple:
    // Minimal stand-in:
    let mut h = key;
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

impl CompactHashTable {
    pub fn new(capacity: usize, key_bits: usize, value_bits: usize) -> Result<Self> {
        if key_bits + value_bits != std::mem::size_of::<CompactHashCell>() * 8 {
            bail!("sum of key bits and value bits must equal 64");
        }
        if key_bits == 0 || value_bits == 0 {
            bail!("key_bits and value_bits cannot be zero");
        }

        let mut zone_locks_ = Vec::with_capacity(LOCK_ZONES);
        for _ in 0..LOCK_ZONES {
            zone_locks_.push(Mutex::new(()));
        }

        let mut backing_data_ = Vec::with_capacity(capacity);
        for _ in 0..capacity {
            backing_data_.push(CompactHashCell::new());
        }

        let table_ptr = backing_data_.as_mut_ptr();

        Ok(CompactHashTable {
            capacity_: capacity,
            size_: 0,
            key_bits_: key_bits,
            value_bits_: value_bits,
            file_backed_: false,
            locks_initialized_: true,
            table_: table_ptr,
            zone_locks_,
            backing_data_,
        })
    }

    pub fn from_file(filename: &str, memory_mapping: bool) -> Result<Self> {
        if memory_mapping {
            // Memory mapping not implemented
            bail!("memory mapping not implemented");
        } else {
            let mut f = File::open(filename)?;
            let mut capacity = 0usize;
            let mut size = 0usize;
            let mut key_bits = 0usize;
            let mut value_bits = 0usize;

            unsafe {
                let buf = std::slice::from_raw_parts_mut(
                    &mut capacity as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                );
                f.read_exact(buf)?;
            }

            unsafe {
                let buf = std::slice::from_raw_parts_mut(
                    &mut size as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                );
                f.read_exact(buf)?;
            }

            unsafe {
                let buf = std::slice::from_raw_parts_mut(
                    &mut key_bits as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                );
                f.read_exact(buf)?;
            }

            unsafe {
                let buf = std::slice::from_raw_parts_mut(
                    &mut value_bits as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                );
                f.read_exact(buf)?;
            }

            let mut backing_data_ = vec![CompactHashCell::new(); capacity];
            unsafe {
                let bytes_needed = capacity * std::mem::size_of::<CompactHashCell>();
                let buf = std::slice::from_raw_parts_mut(
                    backing_data_.as_mut_ptr() as *mut u8,
                    bytes_needed,
                );
                f.read_exact(buf)?;
            }

            let mut zone_locks_ = Vec::with_capacity(LOCK_ZONES);
            for _ in 0..LOCK_ZONES {
                zone_locks_.push(Mutex::new(()));
            }

            let table_ptr = backing_data_.as_mut_ptr();
            Ok(CompactHashTable {
                capacity_: capacity,
                size_: size,
                key_bits_: key_bits,
                value_bits_: value_bits,
                file_backed_: false,
                locks_initialized_: true,
                table_: table_ptr,
                zone_locks_,
                backing_data_,
            })
        }
    }

    pub fn write_table(&self, filename: &str) -> Result<()> {
        let mut ofs = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(filename)?;
        unsafe {
            let buf_cap = std::slice::from_raw_parts(
                &self.capacity_ as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            );
            ofs.write_all(buf_cap)?;

            let buf_size = std::slice::from_raw_parts(
                &self.size_ as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            );
            ofs.write_all(buf_size)?;

            let buf_kb = std::slice::from_raw_parts(
                &self.key_bits_ as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            );
            ofs.write_all(buf_kb)?;

            let buf_vb = std::slice::from_raw_parts(
                &self.value_bits_ as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            );
            ofs.write_all(buf_vb)?;

            let bytes_needed = self.capacity_ * std::mem::size_of::<CompactHashCell>();
            let table_bytes =
                std::slice::from_raw_parts(self.backing_data_.as_ptr() as *const u8, bytes_needed);
            ofs.write_all(table_bytes)?;
        }
        Ok(())
    }

    pub fn get(&self, key: hkey_t) -> hvalue_t {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits_);
        let mut idx = (hc % (self.capacity_ as u64)) as usize;
        let first_idx = idx;
        let mut step = 0;

        // let table = unsafe { std::slice::from_raw_parts(self.table_, self.capacity_) };
        let table = &self.backing_data_; // Use backing_data_ directly instead of unsafe pointers

        loop {
            let val = table[idx].value(self.value_bits_);
            if val == 0 {
                break; // empty cell
            }
            if table[idx].hashed_key(self.value_bits_) == compacted_key {
                return val;
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity_;
            if idx == first_idx {
                break;
            }
        }
        0
    }

    pub fn find_index(&self, key: hkey_t, idx_out: &mut usize) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits_);
        *idx_out = (hc % (self.capacity_ as u64)) as usize;
        let first_idx = *idx_out;
        let mut step = 0;

        let table = unsafe { std::slice::from_raw_parts(self.table_, self.capacity_) };
        loop {
            let val = table[*idx_out].value(self.value_bits_);
            if val == 0 {
                return false; // empty cell
            }
            if table[*idx_out].hashed_key(self.value_bits_) == compacted_key {
                return true;
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            *idx_out = (*idx_out + step) % self.capacity_;
            if *idx_out == first_idx {
                break;
            }
        }
        false
    }

    pub fn compare_and_set(
        &mut self,
        key: hkey_t,
        new_value: hvalue_t,
        old_value: &mut hvalue_t,
    ) -> bool {
        if self.file_backed_ {
            return false;
        }
        if new_value == 0 {
            return false;
        }

        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits_);
        let mut idx = (hc % (self.capacity_ as u64)) as usize;
        let first_idx = idx;
        let mut step = 0;
        let table = unsafe { std::slice::from_raw_parts_mut(self.table_, self.capacity_) };

        loop {
            let zone = idx % LOCK_ZONES;
            let _guard = self.zone_locks_[zone].lock().unwrap();
            let cell_val = table[idx].value(self.value_bits_);
            let cell_key = table[idx].hashed_key(self.value_bits_);
            if cell_val == 0 || cell_key == compacted_key {
                // location found
                if *old_value == cell_val {
                    table[idx].populate(compacted_key, new_value, self.key_bits_, self.value_bits_);
                    if *old_value == 0 {
                        self.size_ += 1;
                    }
                    return true;
                } else {
                    *old_value = cell_val;
                    return false;
                }
            }
            drop(_guard);
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity_;
            if idx == first_idx {
                return false;
            }
        }
    }

    pub fn direct_compare_and_set(
        &mut self,
        idx: usize,
        key: hkey_t,
        new_value: hvalue_t,
        old_value: &mut hvalue_t,
    ) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits_);
        let zone = idx % LOCK_ZONES;
        let table = unsafe { std::slice::from_raw_parts_mut(self.table_, self.capacity_) };
        let _guard = self.zone_locks_[zone].lock().unwrap();
        let cell_val = table[idx].value(self.value_bits_);
        if *old_value == cell_val {
            table[idx].populate(compacted_key, new_value, self.key_bits_, self.value_bits_);
            if *old_value == 0 {
                self.size_ += 1;
            }
            true
        } else {
            *old_value = cell_val;
            false
        }
    }

    fn second_hash(&self, first_hash: u64) -> usize {
        if LINEAR_PROBING {
            1
        } else {
            ((first_hash >> 8) | 1) as usize
        }
    }

    pub fn get_value_counts(&self) -> taxon_counts_t {
        let table = unsafe { std::slice::from_raw_parts(self.table_, self.capacity_) };
        let mut value_counts = HashMap::new();
        for i in 0..self.capacity_ {
            let val = table[i].value(self.value_bits_);
            if val != 0 {
                *value_counts.entry(val).or_insert(0) += 1;
            }
        }
        value_counts
    }

    pub fn size(&self) -> usize {
        self.size_
    }

    pub fn capacity(&self) -> usize {
        self.capacity_
    }
}
