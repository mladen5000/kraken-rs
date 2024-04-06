use crate::hyperloglogplus::HyperLogLogPlusMinus;
use crate::readcounts::ReadCounts;
use std::collections::HashMap;

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

impl IndexOptions {
    pub fn new() -> Self {
        IndexOptions {
            k: 0,
            l: 0,
            spaced_seed_mask: 0,
            toggle_mask: 0,
            dna_db: false,
            minimum_acceptable_hash_value: 0,
            revcom_version: 0,
            db_version: 0,
            db_type: 0,
        }
    }
}

pub type TaxId = u64;
pub const TAXID_MAX: TaxId = !0;

pub type TaxonCounts = HashMap<TaxId, u64>;

#[cfg(feature = "exact_counting")]
pub type ReadCounter = ReadCounts<HashSet<u64>>;
#[cfg(not(feature = "exact_counting"))]
pub type ReadCounter = ReadCounts<HyperLogLogPlusMinus>;

pub type TaxonCounters = HashMap<TaxId, ReadCounter>;
