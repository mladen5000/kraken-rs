/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use serde::Deserialize;
use serde::Serialize;

use crate::readcounts::HyperLogLogPlusMinus;
use crate::readcounts::ReadCounts;
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
        bincode::deserialize(&buffer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
pub type TaxId = u64;
pub const TAXID_MAX: TaxId = u64::MAX;

pub type TaxonCounts = HashMap<TaxId, u64>;

// Conditional compilation based on the `exact_counting` feature
#[cfg(feature = "exact_counting")]
pub type ReadCounter = ReadCounts<HashSet<u64>>;

#[cfg(not(feature = "exact_counting"))]
pub type ReadCounter = ReadCounts<HyperLogLogPlusMinus>;

pub type TaxonCounters = HashMap<TaxId, ReadCounter>;
