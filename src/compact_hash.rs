use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::kv_store::KeyValueStore;

pub struct CompactHashTable {
    table: Box<[AtomicU64]>,
    value_bits: u8,
    taxid_bits: u8,
    value_mask: u64,
    taxid_mask: u64,
}

impl CompactHashTable {
    pub fn new(capacity: usize, value_bits: u8, taxid_bits: u8) -> Self {
        let table = (0..capacity)
            .map(|_| AtomicU64::new(0))
            .collect::<Vec<_>>()
            .into_boxed_slice();

        let value_mask = (1 << value_bits) - 1;
        let taxid_mask = (1 << taxid_bits) - 1;

        Self {
            table,
            value_bits,
            taxid_bits,
            value_mask,
            taxid_mask,
        }
    }

    pub fn write_table(&self, path: &str) -> Result<()> {
        let mut file = File::create(path)?;
        for value in self.table.iter() {
            file.write_all(&value.load(Ordering::Relaxed).to_le_bytes())?;
        }
        Ok(())
    }

    pub fn node_count(&self) -> usize {
        self.table.len()
    }

    pub fn get_at(&self, idx: usize) -> &AtomicU64 {
        &self.table[idx]
    }

    pub fn find_index(&self, key: u64, idx: &mut usize) -> bool {
        let probe_idx = (key & (self.table.len() as u64 - 1)) as usize;
        *idx = probe_idx;
        let value = self.table[probe_idx].load(Ordering::Relaxed);
        let stored_key = (value >> (self.value_bits + self.taxid_bits)) & self.value_mask;
        stored_key == key
    }

    pub fn compare_and_set(&mut self, key: u64, new_value: u64, old_value: &mut u64) -> bool {
        let probe_idx = (key & (self.table.len() as u64 - 1)) as usize;
        let current = self.table[probe_idx].load(Ordering::Relaxed);
        let stored_key = (current >> (self.value_bits + self.taxid_bits)) & self.value_mask;

        if stored_key == 0 {
            // Empty slot - we can insert
            *old_value = 0;
            let new_packed =
                (key << (self.value_bits + self.taxid_bits)) | (new_value & self.taxid_mask);

            self.table[probe_idx].store(new_packed, Ordering::Relaxed);
            true
        } else if stored_key == key {
            // Found existing key
            *old_value = current & self.taxid_mask;
            false
        } else {
            // Different key in this slot
            *old_value = 0;
            false
        }
    }

    pub fn set_minimizer_lca(
        &mut self,
        key: u64,
        taxid: u64,
        taxonomy: &crate::taxonomy::Taxonomy,
    ) {
        let mut old_value = 0;
        let mut new_value = taxid;
        while !self.compare_and_set(key, new_value, &mut old_value) {
            if old_value == 0 {
                // Means there was a different key in the slot
                break;
            }
            new_value = taxonomy.lowest_common_ancestor(old_value, taxid);
        }
    }

    pub fn size(&self) -> usize {
        let mut count = 0;
        for value in self.table.iter() {
            if value.load(Ordering::Relaxed) != 0 {
                count += 1;
            }
        }
        count
    }

    pub fn capacity(&self) -> usize {
        self.table.len()
    }

    pub fn get_value_counts(&self) -> HashMap<u64, u64> {
        let mut counts = HashMap::new();
        for value in self.table.iter() {
            let packed = value.load(Ordering::Relaxed);
            if packed != 0 {
                let taxid = packed & self.taxid_mask;
                *counts.entry(taxid).or_insert(0) += 1;
            }
        }
        counts
    }
}

impl KeyValueStore for CompactHashTable {
    fn get(&self, key: u64) -> u64 {
        let mut idx = 0;
        if self.find_index(key, &mut idx) {
            self.table[idx].load(Ordering::Relaxed) & self.taxid_mask
        } else {
            0
        }
    }
}
