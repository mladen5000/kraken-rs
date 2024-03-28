// pub mod aa_translate;
// pub mod build_db;
// pub mod classify;
// pub mod classify3;
// pub mod compact_hash;
// mod estimate_capacity;
// pub mod hyperloglogplus;
// pub mod kraken2_data;
// pub mod kv_store;
// pub mod mmap_file;
// pub mod mmscanner;
// pub mod readcounts;
// pub mod reports;
pub mod seqreader;
// pub mod taxonomy;
// pub mod threadpool;
// pub mod utilities;
pub mod k2mask;

use k2mask::main as k2mask_main;

fn main() {
    println!("Hello, world!");
    k2mask_main();
}
