use crate::hyperloglogplus::HyperLogLogPlusMinus;
use crate::reports::ReadCounts;
use std::collections::HashMap;
use std::collections::HashSet;

pub struct IndexOptions {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    dna_db: bool,
    pub minimum_acceptable_hash_value: u64,
    revcom_version: i32, // Fix bug from before K2.0.8
    db_version: i32,     // To allow for future database structural changes
    db_type: i32,        // To allow for future use of other data structures
}

pub type TaxId = u64;
pub const TAXID_MAX: TaxId = !0;

pub type TaxonCounts = HashMap<TaxId, u64>;

#[cfg(feature = "exact_counting")]
pub type ReadCounter = ReadCounts<HashSet<u64>>;
#[cfg(not(feature = "exact_counting"))]
pub type ReadCounter = ReadCounts<HyperLogLogPlusMinus<u64>>;

pub type TaxonCounters = HashMap<TaxId, ReadCounter>;
