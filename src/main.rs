pub mod aa_translate;
pub mod classify;
pub mod estimate_capacity;
pub mod hyperloglogplus;

pub mod kraken2_data;
pub mod kv_store;
pub mod mladen;
pub mod reports;
pub mod threadpool;

// Build DB
pub mod build_db;
pub mod compact_hash;
pub mod kraken2_headers;
pub mod mmscanner;
pub mod omp_hack;
pub mod seqreader;
pub mod taxonomy;
pub mod utilities;

//  lookup accession numbers
pub mod lookup_accession_numbers;
pub mod readcounts;

// k2mask
pub mod gz_stream;
pub mod k2mask;
fn main() {
    println!("Hello, world!");

    k2mask::run();
}
