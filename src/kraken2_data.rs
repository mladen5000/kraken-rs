/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use crate::readcounts::HyperLogLogPlusMinus;
use crate::readcounts::ReadCounts;
use std::collections::HashMap;
use std::collections::HashSet;

pub struct IndexOptions {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub dna_db: bool,
    pub minimum_acceptable_hash_value: u64,
    pub revcom_version: i32, // Fix bug from before K2.0.8
    pub db_version: i32,     // To allow for future database structural changes
    pub db_type: i32,        // To allow for future use of other data structures
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
