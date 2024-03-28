/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use crate::kraken2_data::TaxonCounts;
use crate::kv_store::{murmur_hash3, HKey, HValue, KeyValueStore};
use crate::mmap_file::MMapFile;

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use std::fs::File;
use std::io::{Read, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::{mem, ptr};

const LOCK_ZONES: usize = 256;

#[derive(Debug, Default)]
struct CompactHashCell {
    data: u32,
}

impl CompactHashCell {
    fn hashed_key(&self, value_bits: usize) -> HKey {
        (self.data >> value_bits) as HKey
    }

    fn value(&self, value_bits: usize) -> HValue {
        (self.data & ((1 << value_bits) - 1)) as HValue
    }

    fn populate(&mut self, compacted_key: HKey, val: HValue, key_bits: usize, value_bits: usize) {
        if key_bits + value_bits != 32 {
            panic!(
                "key len of {} and value len of {} don't sum to 32",
                key_bits, value_bits
            );
        }
        if key_bits == 0 || value_bits == 0 {
            panic!("key len and value len must be nonzero");
        }
        let max_value = (1u64 << value_bits) - 1;
        if max_value < val as u64 {
            panic!("value len of {} too small for value of {}", value_bits, val);
        }
        self.data = ((compacted_key as u32) << value_bits) | (val as u32);
    }
}

pub struct CompactHashTable {
    capacity: usize,
    size: AtomicUsize,
    key_bits: usize,
    value_bits: usize,
    table: Vec<CompactHashCell>,
    file_backed: bool,
    locks_initialized: bool,
    backing_file: Option<MMapFile>,
    zone_locks: Vec<Mutex<()>>,
}

impl CompactHashTable {
    pub fn new(capacity: usize, key_bits: usize, value_bits: usize) -> Self {
        if key_bits + value_bits != mem::size_of::<CompactHashCell>() * 8 {
            panic!(
                "sum of key bits and value bits must equal {}",
                mem::size_of::<CompactHashCell>() * 8
            );
        }
        if key_bits == 0 {
            panic!("key bits cannot be zero");
        }
        if value_bits == 0 {
            panic!("value bits cannot be zero");
        }
        let mut table = Vec::with_capacity(capacity);
        table.resize_with(capacity, CompactHashCell::default);
        let zone_locks = std::iter::repeat_with(|| Mutex::new(()))
            .take(LOCK_ZONES)
            .collect();
        Self {
            capacity,
            size: AtomicUsize::new(0),
            key_bits,
            value_bits,
            table,
            file_backed: false,
            locks_initialized: true,
            backing_file: None,
            zone_locks,
        }
    }

    pub fn from_file(filename: &str, memory_mapping: bool) -> Self {
        let mut cht = Self::default();
        cht.load_table(filename, memory_mapping);
        cht
    }

    fn default() -> Self {
        Self {
            capacity: 0,
            size: AtomicUsize::new(0),
            key_bits: 0,
            value_bits: 0,
            table: Vec::new(),
            file_backed: false,
            locks_initialized: false,
            backing_file: None,
            zone_locks: Vec::new(),
        }
    }

    fn load_table(&mut self, filename: &str, memory_mapping: bool) {
        self.locks_initialized = false;
        if memory_mapping {
            let mut backing_file = MMapFile::new();
            backing_file.open_file(filename, usize::MAX);
            let ptr = backing_file.fptr();
            let mut offset = 0;
            self.capacity = unsafe { *(ptr.add(offset) as *const usize) };
            offset += mem::size_of::<usize>();
            self.size = AtomicUsize::new(unsafe { *(ptr.add(offset) as *const usize) });
            offset += mem::size_of::<usize>();
            self.key_bits = unsafe { *(ptr.add(offset) as *const usize) };
            offset += mem::size_of::<usize>();
            self.value_bits = unsafe { *(ptr.add(offset) as *const usize) };
            offset += mem::size_of::<usize>();
            self.table = unsafe {
                let slice = std::slice::from_raw_parts(
                    ptr.add(offset) as *const CompactHashCell,
                    self.capacity,
                );
                (0..self.capacity).map(|i| ptr::read(&slice[i])).collect()
            };
            if backing_file.filesize() - offset != mem::size_of::<CompactHashCell>() * self.capacity
            {
                panic!("Capacity mismatch in {}, aborting", filename);
            }
            self.file_backed = true;
            self.backing_file = Some(backing_file);
        } else {
            let mut file = File::open(filename).unwrap();
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut self.capacity as *mut usize as *mut u8,
                    mem::size_of::<usize>(),
                )
            })
            .unwrap();
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut self.size as *mut AtomicUsize as *mut u8,
                    mem::size_of::<usize>(),
                )
            })
            .unwrap();
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut self.key_bits as *mut usize as *mut u8,
                    mem::size_of::<usize>(),
                )
            })
            .unwrap();
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut self.value_bits as *mut usize as *mut u8,
                    mem::size_of::<usize>(),
                )
            })
            .unwrap();
            self.table = Vec::with_capacity(self.capacity);
            for _ in 0..self.capacity {
                self.table.push(CompactHashCell { data: 0 });
            }
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    self.table.as_mut_ptr() as *mut u8,
                    mem::size_of::<CompactHashCell>() * self.capacity,
                )
            })
            .unwrap();
            self.file_backed = false;
        }
    }

    pub fn write_table(&self, filename: &str) {
        let mut file = File::create(filename).unwrap();
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.capacity as *const usize as *const u8,
                mem::size_of::<usize>(),
            )
        })
        .unwrap();
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.size as *const AtomicUsize as *const u8,
                mem::size_of::<usize>(),
            )
        })
        .unwrap();
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.key_bits as *const usize as *const u8,
                mem::size_of::<usize>(),
            )
        })
        .unwrap();
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.value_bits as *const usize as *const u8,
                mem::size_of::<usize>(),
            )
        })
        .unwrap();
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                self.table.as_ptr() as *const u8,
                mem::size_of::<CompactHashCell>() * self.capacity,
            )
        })
        .unwrap();
    }

    fn second_hash(first_hash: u64) -> u64 {
        // Linear probing
        1
    }

    pub fn get_value_counts(&self) -> TaxonCounts {
        let mut value_counts = TaxonCounts::new();
        let num_threads = num_cpus::get();
        let mut thread_value_counts = vec![TaxonCounts::new(); num_threads];
        let value_bits = self.value_bits;

        rayon::iter::IndexedParallelIterator::enumerate(self.table.par_iter()).for_each(
            |(i, cell)| {
                let val = cell.value(value_bits);
                if val != 0 {
                    thread_value_counts[i % num_threads]
                        .entry(val.into())
                        .and_modify(|count| *count += 1)
                        .or_insert(1);
                }
            },
        );
        for thread_counts in thread_value_counts {
            for (val, count) in thread_counts {
                value_counts
                    .entry(val)
                    .and_modify(|c| *c += count)
                    .or_insert(count);
            }
        }
        value_counts
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn size(&self) -> usize {
        self.size.load(Ordering::SeqCst)
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

    fn find_index(&self, key: HKey, idx: &mut usize) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        *idx = (hc % self.capacity as u64) as usize;
        let first_idx = *idx;
        let mut step = 0;
        loop {
            if self.table[*idx].value(self.value_bits) == 0 {
                return false;
            }
            if self.table[*idx].hashed_key(self.value_bits) == compacted_key {
                return true;
            }
            if step == 0 {
                step = Self::second_hash(hc);
            }
            *idx += step as usize;
            *idx %= self.capacity;
            if *idx == first_idx {
                break;
            }
        }
        false
    }

    fn compare_and_set(&self, key: HKey, new_value: HValue, old_value: &mut HValue) -> bool {
        if self.file_backed {
            return false;
        }
        if new_value == 0 {
            return false;
        }
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step = 0;
        let mut set_successful = false;
        let mut search_successful = false;
        while !search_successful {
            let zone = idx % LOCK_ZONES;
            let _lock = self.zone_locks[zone].lock().unwrap();
            if self.table[idx].value(self.value_bits) == 0
                || self.table[idx].hashed_key(self.value_bits) == compacted_key
            {
                search_successful = true;
                if *old_value == self.table[idx].value(self.value_bits) {
                    self.table[idx].populate(
                        compacted_key,
                        new_value,
                        self.key_bits,
                        self.value_bits,
                    );
                    if *old_value == 0 {
                        self.size.fetch_add(1, Ordering::SeqCst);
                    }
                    set_successful = true;
                } else {
                    *old_value = self.table[idx].value(self.value_bits);
                }
            }
            if step == 0 {
                step = Self::second_hash(hc);
            }
            idx += step as usize;
            idx %= self.capacity;
            if idx == first_idx {
                panic!("compact hash table capacity exceeded");
            }
        }
        set_successful
    }

    fn direct_compare_and_set(
        &self,
        idx: usize,
        key: HKey,
        new_value: HValue,
        old_value: &mut HValue,
    ) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut set_successful = false;
        let zone = idx % LOCK_ZONES;
        let _lock = self.zone_locks[zone].lock().unwrap();
        if *old_value == self.table[idx].value(self.value_bits) {
            self.table[idx].populate(compacted_key, new_value, self.key_bits, self.value_bits);
            if *old_value == 0 {
                self.size.fetch_add(1, Ordering::SeqCst);
            }
            set_successful = true;
        } else {
            *old_value = self.table[idx].value(self.value_bits);
        }
        set_successful
    }
}

impl KeyValueStore for CompactHashTable {
    fn get(&self, key: HKey) -> HValue {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = (hc % self.capacity as u64) as usize;
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
                step = Self::second_hash(hc);
            }
            idx += step as usize;
            idx %= self.capacity;
            if idx == first_idx {
                break;
            }
        }
        0
    }
}
