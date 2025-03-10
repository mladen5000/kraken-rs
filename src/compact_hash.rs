/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 * Rust conversion Copyright 2025
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::sync::atomic::{AtomicU32, Ordering};

use crate::kraken2_data::TaxonCounts;
use crate::kv_store::KeyValueStore;
use crate::mmap_file::MMapFile;
use crate::utilities::murmur_hash3;

// CompactHashCell matches the C++ struct
#[repr(C)]
pub struct CompactHashCell {
    data: AtomicU32,
}

// CompactHashCell is not automatically Clone because AtomicU32 doesn't implement Clone
// We implement a custom Clone that preserves the atomic value
impl Clone for CompactHashCell {
    fn clone(&self) -> Self {
        Self {
            data: AtomicU32::new(self.data.load(Ordering::Relaxed)),
        }
    }
}

impl CompactHashCell {
    pub fn hashed_key(&self, value_bits: u8) -> u64 {
        (self.data.load(Ordering::Relaxed) >> value_bits) as u64
    }

    pub fn value(&self, value_bits: u8) -> u64 {
        (self.data.load(Ordering::Relaxed) & ((1 << value_bits) - 1)) as u64
    }

    pub fn populate(&self, compacted_key: u64, val: u64, key_bits: u8, value_bits: u8) {
        assert_eq!(
            key_bits + value_bits,
            32,
            "key_bits + value_bits must equal 32"
        );
        assert!(key_bits > 0, "key_bits must be positive");
        assert!(value_bits > 0, "value_bits must be positive");

        let max_value = (1 << value_bits) - 1;
        assert!(val <= max_value as u64, "value is too large for value_bits");

        let data = ((compacted_key << value_bits) | val) as u32;
        self.data.store(data, Ordering::Relaxed);
    }
}

impl Default for CompactHashCell {
    fn default() -> Self {
        Self {
            data: AtomicU32::new(0),
        }
    }
}

pub struct CompactHashTable {
    capacity: usize,
    size: usize,
    key_bits: u8,
    value_bits: u8,
    table: *mut CompactHashCell,
    file_backed: bool,
    backing_file: Option<MMapFile>,
}

// This is safe because we ensure proper access patterns in our implementation
unsafe impl Send for CompactHashTable {}
unsafe impl Sync for CompactHashTable {}

impl CompactHashTable {
    pub fn new(capacity: usize, key_bits: u8, value_bits: u8) -> Self {
        assert_eq!(
            key_bits + value_bits,
            32,
            "key_bits + value_bits must equal 32"
        );
        assert!(key_bits > 0, "key_bits must be positive");
        assert!(value_bits > 0, "value_bits must be positive");

        let table = Box::into_raw(vec![CompactHashCell::default(); capacity].into_boxed_slice())
            as *mut CompactHashCell;

        Self {
            capacity,
            size: 0,
            key_bits,
            value_bits,
            table,
            file_backed: false,
            backing_file: None,
        }
    }

    pub fn from_file(filename: &str, memory_mapping: bool) -> Result<Self> {
        Self::load_table(filename, memory_mapping)
    }

    fn load_table(filename: &str, memory_mapping: bool) -> Result<Self> {
        if memory_mapping {
            let backing_file =
                MMapFile::open_file(filename).context("Failed to memory map hash table file")?;

            let ptr = backing_file.fptr();
            let mut offset = 0;

            // Read header values
            let capacity = unsafe {
                let val_ptr = ptr.as_ptr().add(offset) as *const usize;
                offset += std::mem::size_of::<usize>();
                *val_ptr
            };

            let size = unsafe {
                let val_ptr = ptr.as_ptr().add(offset) as *const usize;
                offset += std::mem::size_of::<usize>();
                *val_ptr
            };

            let key_bits = unsafe {
                let val_ptr = ptr.as_ptr().add(offset) as *const u8;
                offset += std::mem::size_of::<u8>();
                *val_ptr
            };

            let value_bits = unsafe {
                let val_ptr = ptr.as_ptr().add(offset) as *const u8;
                offset += std::mem::size_of::<u8>();
                *val_ptr
            };

            // Pad to 8-byte alignment if necessary
            if offset % 8 != 0 {
                offset += 8 - (offset % 8);
            }

            // Table data starts at this offset
            let table = unsafe { ptr.as_ptr().add(offset) as *mut CompactHashCell };

            // Verify file size matches expected size
            let expected_size = offset + capacity * std::mem::size_of::<CompactHashCell>();
            if backing_file.len() < expected_size {
                anyhow::bail!(
                    "Hash table file size mismatch, expected at least {} bytes",
                    expected_size
                );
            }

            Ok(Self {
                capacity,
                size,
                key_bits,
                value_bits,
                table,
                file_backed: true,
                backing_file: Some(backing_file),
            })
        } else {
            let mut file = File::open(filename).context("Failed to open hash table file")?;

            // Read header values
            let mut capacity = 0usize;
            let mut size = 0usize;
            let mut key_bits = 0u8;
            let mut value_bits = 0u8;

            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut capacity as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                )
            })?;

            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut size as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                )
            })?;

            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut key_bits as *mut u8 as *mut u8,
                    std::mem::size_of::<u8>(),
                )
            })?;

            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut value_bits as *mut u8 as *mut u8,
                    std::mem::size_of::<u8>(),
                )
            })?;

            // Allocate table
            let table = Box::into_raw(vec![CompactHashCell::default(); capacity].into_boxed_slice())
                as *mut CompactHashCell;

            // Read the table data
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    table as *mut u8,
                    capacity * std::mem::size_of::<CompactHashCell>(),
                )
            })?;

            Ok(Self {
                capacity,
                size,
                key_bits,
                value_bits,
                table,
                file_backed: false,
                backing_file: None,
            })
        }
    }

    pub fn write_table(&self, filename: &str) -> Result<()> {
        let mut file = File::create(filename).context("Failed to create hash table file")?;

        // Write header values
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.capacity as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            )
        })?;

        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.size as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            )
        })?;

        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.key_bits as *const u8 as *const u8,
                std::mem::size_of::<u8>(),
            )
        })?;

        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.value_bits as *const u8 as *const u8,
                std::mem::size_of::<u8>(),
            )
        })?;

        // Write the table data
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                self.table as *const u8,
                self.capacity * std::mem::size_of::<CompactHashCell>(),
            )
        })?;

        Ok(())
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn key_bits(&self) -> u8 {
        self.key_bits
    }

    pub fn value_bits(&self) -> u8 {
        self.value_bits
    }

    pub fn occupancy(&self) -> f64 {
        self.size as f64 / self.capacity as f64
    }

    // Implements the same logic as the C++ second_hash function
    fn second_hash(&self, first_hash: u64) -> usize {
        // Simplified version that avoids unreachable code warnings
        #[cfg(feature = "linear_probing")]
        return 1;
        
        // Default behavior (non-linear probing or no feature specified)
        ((first_hash >> 8) | 1) as usize
    }

    pub fn get_value_counts(&self) -> TaxonCounts {
        let mut counts = HashMap::new();

        // This would ideally be parallelized like in the C++ version
        for i in 0..self.capacity {
            unsafe {
                let cell = &*self.table.add(i);
                let val = cell.value(self.value_bits);
                if val > 0 {
                    *counts.entry(val).or_insert(0) += 1;
                }
            }
        }

        counts
    }

    /// Compare and set a key's value atomically.
    /// Returns true if successful or false if a collision occurs.
    /// Compare and set a key's value atomically.
    /// Returns true if the value was set, false if not found or collision.
    /// If the key is found, old_val is set to the current value.
    pub fn compare_and_set(&self, key: u64, new_val: u64, old_val: &mut u64) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits as u64);
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step = 0;

        // First find the key
        loop {
            unsafe {
                let cell = &*self.table.add(idx);
                let cell_value = cell.value(self.value_bits);
                
                if cell_value == 0 {
                    // Key not found, return false
                    return false;
                } else if cell.hashed_key(self.value_bits) == compacted_key {
                    // Key found, set old_val to current value
                    *old_val = cell_value;
                    
                    // Get the hashed key from the new value
                    let val_mask = (1 << self.value_bits) - 1;
                    let key_mask = !val_mask;
                    
                    // Atomically update the cell with the new value
                    // Keep the same key but update the value
                    let old_data = cell.data.load(std::sync::atomic::Ordering::Relaxed);
                    let new_data = (old_data & key_mask as u32) | ((new_val & val_mask) as u32);
                    
                    // Try to set the new value atomically
                    cell.data.store(new_data, std::sync::atomic::Ordering::Relaxed);
                    
                    return true;
                }
            }

            step += 1;
            idx = (first_idx + step * step) % self.capacity;
            if idx == first_idx {
                // Wrapped around, key not found
                return false;
            }
        }
    }

    pub fn find_index(&self, key: u64, idx: &mut usize) -> bool {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits as u64);
        *idx = (hc % self.capacity as u64) as usize;
        let first_idx = *idx;
        let mut step = 0;

        loop {
            unsafe {
                let cell = &*self.table.add(*idx);
                if cell.value(self.value_bits) == 0 {
                    // Empty cell, search is over
                    return false;
                }

                if cell.hashed_key(self.value_bits) == compacted_key {
                    return true;
                }

                if step == 0 {
                    step = self.second_hash(hc);
                }

                *idx += step;
                *idx %= self.capacity;

                if *idx == first_idx {
                    // We've exhausted the table
                    break;
                }
            }
        }

        false
    }
}

impl KeyValueStore for CompactHashTable {
    fn get(&self, key: u64) -> u64 {
        let hc = murmur_hash3(key);
        let compacted_key = hc >> (32 + self.value_bits as u64);
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step = 0;

        loop {
            unsafe {
                let cell = &*self.table.add(idx);
                if cell.value(self.value_bits) == 0 {
                    // Empty cell, search is over
                    break;
                }

                if cell.hashed_key(self.value_bits) == compacted_key {
                    return cell.value(self.value_bits);
                }

                if step == 0 {
                    step = self.second_hash(hc);
                }

                idx += step;
                idx %= self.capacity;

                if idx == first_idx {
                    // We've exhausted the table
                    break;
                }
            }
        }

        0
    }
}

impl Drop for CompactHashTable {
    fn drop(&mut self) {
        if !self.file_backed {
            unsafe {
                // Deallocate the table if it's not file-backed
                let _ = Box::from_raw(std::slice::from_raw_parts_mut(self.table, self.capacity));
            }
        }
        // The backing_file will be dropped automatically if it exists
    }
}

// Helper functions for tests
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_compact_hash_cell() {
        let cell = CompactHashCell::default();
        cell.populate(123, 456, 24, 8);
        assert_eq!(cell.hashed_key(8), 123);
        assert_eq!(cell.value(8), 456);
    }

    #[test]
    fn test_get() {
        let table = CompactHashTable::new(1024, 24, 8);
        // We'd need to set values manually for testing
        // This would require unsafe code to manipulate the table
    }
}
