// build_db and classify are now in bin/ as separate binaries
pub mod compact_hash;
pub mod mmap_file;
pub mod mmscanner;
pub mod omp_hack;
pub mod seqreader;
pub mod taxonomy;
pub mod utilities;

// classify
pub mod aa_translate;
pub mod hyperloglogplus;
pub mod reports;

//estimate_capacity

// dump_table

// k2mask
pub mod kraken2_data;
pub mod kv_store;
pub mod threadpool;

//  lookup session numbers
pub mod readcounts;

// k2mask
pub mod gz_stream;
/// Main entry point for the k2mask executable.
///
/// Prints "Hello, world!" to STDOUT, then runs the k2mask::run() function.
fn main() {
    println!("Hello, world!");
    bin::build_db::main();
}
