/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use serde::Deserialize;
use serde::Serialize;

use crate::hyperloglogplus::HyperLogLogPlusMinus;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Read;
use std::path::Path;

#[derive(Default, Serialize, Deserialize)]
pub struct IndexOptions {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub dna_db: bool,
    pub minimum_acceptable_hash_value: u64,
    pub revcom_version: u32, // Fix bug from before K2.0.8
    pub db_version: u32,     // To allow for future database structural changes
    pub db_type: u32,        // To allow for future use of other data structures
}

impl IndexOptions {
    pub fn load(filename: &Path) -> io::Result<Self> {
        let mut file = File::open(filename)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        // Use our own implementation since we just need the bytes copied
        let size = std::mem::size_of::<Self>();
        if buffer.len() < size {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "buffer too small"));
        }
        
        let mut result = Self::default();
        unsafe {
            let dst_ptr = &mut result as *mut Self as *mut u8;
            let src_ptr = buffer.as_ptr();
            std::ptr::copy_nonoverlapping(src_ptr, dst_ptr, size);
        }
        Ok(result)
    }
}
pub type TaxId = u64;
pub const TAXID_MAX: TaxId = u64::MAX;

// Unified TaxonCounters struct that matches C++ implementation
#[derive(Clone)]
pub struct TaxonCounters {
    pub read_count: u64,
    pub kmer_count: HyperLogLogPlusMinus,
}

impl TaxonCounters {
    pub fn new(read_count: u64, kmer_count: HyperLogLogPlusMinus) -> Self {
        Self {
            read_count,
            kmer_count,
        }
    }

    pub fn new_with_precision(precision: u8) -> Self {
        let kmer_count = HyperLogLogPlusMinus::new(precision, true, crate::utilities::murmur_hash3);
        Self {
            read_count: 0,
            kmer_count,
        }
    }

    pub fn increment_read_count(&mut self) {
        self.read_count += 1;
    }

    pub fn get_read_count(&self) -> u64 {
        self.read_count
    }

    pub fn get_kmer_distinct(&self) -> f64 {
        self.kmer_count.cardinality() as f64
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.kmer_count.insert(kmer);
    }

    pub fn merge(&mut self, other: &TaxonCounters) {
        self.read_count += other.read_count;
        self.kmer_count.merge_copy(&other.kmer_count);
    }
}

impl std::ops::AddAssign for TaxonCounters {
    fn add_assign(&mut self, other: Self) {
        self.read_count += other.read_count;
        self.kmer_count.merge_move(other.kmer_count);
    }
}

pub type TaxonCounts = HashMap<TaxId, u64>;
pub type TaxonCountersMap = HashMap<TaxId, TaxonCounters>;
pub type TaxonCountsMap = HashMap<TaxId, u32>;

pub const BITS_PER_CHAR_DNA: u8 = 2;
pub const BITS_PER_CHAR_PRO: u8 = 5;
