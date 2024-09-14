use memmap2::MmapOptions;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex};

const LOCK_ZONES: usize = 256;

pub struct CompactHashCell {
    data: u32,
}

impl CompactHashCell {
    pub fn hashed_key(&self, value_bits: usize) -> u64 {
        (self.data >> value_bits) as u64
    }

    pub fn value(&self, value_bits: usize) -> u64 {
        (self.data & ((1 << value_bits) - 1)) as u64
    }

    pub fn populate(&mut self, compacted_key: u64, val: u64, key_bits: usize, value_bits: usize) {
        if key_bits + value_bits != 32 {
            panic!(
                "key len of {} and value len of {} don't sum to 32",
                key_bits, value_bits
            );
        }
        if key_bits == 0 || value_bits == 0 {
            panic!("key len and value len must be nonzero");
        }
        let max_value = (1 << value_bits) - 1;
        if max_value < val {
            panic!("value len of {} too small for value of {}", value_bits, val);
        }
        self.data = (compacted_key << value_bits) as u32;
        self.data |= val as u32;
    }
}

pub struct CompactHashTable {
    capacity: usize,
    size: Arc<Mutex<usize>>,
    key_bits: usize,
    value_bits: usize,
    table: Vec<CompactHashCell>,
    file_backed: bool,
    locks_initialized: bool,
    zone_locks: Vec<Mutex<()>>,
    backing_file: Option<MmapFile>,
}

impl CompactHashTable {
    pub fn new(capacity: usize, key_bits: usize, value_bits: usize) -> io::Result<Self> {
        if key_bits + value_bits != std::mem::size_of::<CompactHashCell>() * 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "sum of key bits and value bits must equal {}",
                    std::mem::size_of::<CompactHashCell>() * 8
                ),
            ));
        }
        if key_bits == 0 || value_bits == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "key bits and value bits cannot be zero",
            ));
        }

        let mut zone_locks = Vec::with_capacity(LOCK_ZONES);
        for _ in 0..LOCK_ZONES {
            zone_locks.push(Mutex::new(()));
        }

        let table = vec![CompactHashCell { data: 0 }; capacity];

        Ok(CompactHashTable {
            capacity,
            size: Arc::new(Mutex::new(0)),
            key_bits,
            value_bits,
            table,
            file_backed: false,
            locks_initialized: true,
            zone_locks,
            backing_file: None,
        })
    }

    pub fn load_table(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        let mut table = CompactHashTable {
            capacity: 0,
            size: Arc::new(Mutex::new(0)),
            key_bits: 0,
            value_bits: 0,
            table: Vec::new(),
            file_backed: memory_mapping,
            locks_initialized: false,
            zone_locks: vec![Mutex::new(()); LOCK_ZONES],
            backing_file: None,
        };

        if memory_mapping {
            let file = File::open(filename)?;
            let mmap = unsafe { MmapOptions::new().map(&file)? };
            let mut ptr = 0;

            table.capacity = u64::from_ne_bytes(mmap[ptr..ptr + 8].try_into().unwrap()) as usize;
            ptr += 8;
            *table.size.lock().unwrap() = u64::from_ne_bytes(mmap[ptr..ptr + 8].try_into().unwrap()) as usize;
            ptr += 8;
            table.key_bits = u64::from_ne_bytes(mmap[ptr..ptr + 8].try_into().unwrap()) as usize;
            ptr += 8;
            table.value_bits = u64::from_ne_bytes(mmap[ptr..ptr + 8].try_into().unwrap()) as usize;
            ptr += 8;

            let table_size = std::mem::size_of::<CompactHashCell>() * table.capacity;
            if mmap.len() - ptr != table_size {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Capacity mismatch in file",
                ));
            }

            table.table = unsafe {
                std::slice::from_raw_parts(
                    mmap[ptr..].as_ptr() as *const CompactHashCell,
                    table.capacity,
                )
            }
            .to_vec();

            table.backing_file = Some(MmapFile { mmap });
        } else {
            let mut file = File::open(filename)?;
            let mut buffer = [0u8; 8];

            file.read_exact(&mut buffer)?;
            table.capacity = u64::from_ne_bytes(buffer) as usize;

            file.read_exact(&mut buffer)?;
            *table.size.lock().unwrap() = u64::from_ne_bytes(buffer) as usize;

            file.read_exact(&mut buffer)?;
            table.key_bits = u64::from_ne_bytes(buffer) as usize;

            file.read_exact(&mut buffer)?;
            table.value_bits = u64::from_ne_bytes(buffer) as usize;

            table.table = vec![CompactHashCell { data: 0 }; table.capacity];
            for cell in table.table.iter_mut() {
                file.read_exact(unsafe {
                    std::slice::from_raw_parts_mut(
                        &mut cell.data as *mut u32 as *mut u8,
                        std::mem::size_of::<CompactHashCell>(),
                    )
                })?;
            }
        }

        Ok(table)
    }

    pub fn write_table(&self, filename: &str) -> io::Result<()> {
        let mut file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(filename)?;

        file.write_all(&(self.capacity as u64).to_ne_bytes())?;
        file.write_all(&(*self.size.lock().unwrap() as u64).to_ne_bytes())?;
        file.write_all(&(self.key_bits as u64).to_ne_bytes())?;
        file.write_all(&(self.value_bits as u64).to_ne_bytes())?;

        for cell in &self.table {
            file.write_all(&cell.data.to_ne_bytes())?;
        }

        Ok(())
    }

    pub fn get(&self, key: u64) -> u64 {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity;
        let first_idx = idx;
        let mut step = 0;

        loop {
            if self.table[idx].value(self.value_bits) == 0 {
                break;
            }
            if self.table[idx].hashed_key(self.value_bits) == compacted_key {
                return self.table[idx].value(self.value_bits);
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }

        0
    }

    pub fn find_index(&self, key: u64) -> Option<usize> {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity;
        let first_idx = idx;
        let mut step = 0;

        loop {
            if self.table[idx].value(self.value_bits) == 0 {
                return None;
            }
            if self.table[idx].hashed_key(self.value_bits) == compacted_key {
                return Some(idx);
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }

        None
    }

    pub fn compare_and_set(&self, key: u64, new_value: u64, old_value: &mut u64) -> bool {
        if self.file_backed || new_value == 0 {
            return false;
        }
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity;
        let first_idx = idx;
        let mut step = 0;

        loop {
            let zone = idx % LOCK_ZONES;
            let _lock = self.zone_locks[zone].lock().unwrap();

            if self.table[idx].value(self.value_bits) == 0
                || self.table[idx].hashed_key(self.value_bits) == compacted_key
            {
                if *old_value == self.table[idx].value(self.value_bits) {
                    self.table[idx].populate(
                        compacted_key,
                        new_value,
                        self.key_bits,
                        self.value_bits,
                    );
                    if *old_value == 0 {
                        *self.size.lock().unwrap() += 1;
                    }
                    return true;
                } else {
                    *old_value = self.table[idx].value(self.value_bits);
                }
            }

            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                panic!("Compact hash table capacity exceeded");
            }
        }
    }

    pub fn direct_compare_and_set(
        &self,
        idx: usize,
        key: u64,
        new_value: u64,
        old_value: &mut u64,
    ) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let zone = idx % LOCK_ZONES;
        let _lock = self.zone_locks[zone].lock().unwrap();

        if *old_value == self.table[idx].value(self.value_bits) {
            self.table[idx].populate(compacted_key, new_value, self.key_bits, self.value_bits);
            if *old_value == 0 {
                *self.size.lock().unwrap() += 1;
            }
            true
        } else {
            *old_value = self.table[idx].value(self.value_bits);
            false
        }
    }

    pub fn get_value_counts(&self) -> HashMap<u64, usize> {
        let thread_ct = rayon::current_num_threads();
        let chunk_size = (self.capacity + thread_ct - 1) / thread_ct;

        let thread_value_counts: Vec<HashMap<u64, usize>> = self
            .table
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut counts = HashMap::new();
                for cell in chunk {
                    let val = cell.value(self.value_bits);
                    if val != 0 {
                        *counts.entry(val).or_insert(0) += 1;
                    }
                }
                counts
            })
            .collect();

        let mut value_counts = HashMap::new();
        for thread_counts in thread_value_counts {
            for (key, count) in thread_counts {
                *value_counts.entry(key).or_insert(0) += count;
            }
        }
        value_counts
    }

    fn second_hash(&self, first_hash: u64) -> u64 {
        #[cfg(feature = "linear_probing")]
        {
            1
        }
        #[cfg(not(feature = "linear_probing"))]
        {
            (first_hash >> 8) | 1
        }
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn size(&self) -> usize {
        *self.size.lock().unwrap()
    }

    pub fn key_bits(&self) -> usize {
        self.key_bits
    }

    pub fn value_bits(&self) -> usize {
        self.value_bits
    }

    pub fn occupancy(&self) -> f64 {
        self.size() as f64 / self.capacity as f64
    }
}

pub struct MmapFile {
    mmap: memmap2::Mmap,
}

fn murmur_hash3(key: u64) -> u64 {
    // Implement MurmurHash3 here
    // This is a placeholder implementation
    let mut h = key;
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}
