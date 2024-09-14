// pub mod aa_translate;
// pub mod classify;
// pub mod classify3;
// pub mod estimate_capacity;
// pub mod hyperloglogplus;

// pub mod kraken2_data;
// pub mod kv_store;
// // pub mod readcounts_old;
// pub mod reports;
// pub mod threadpool;

// // Build DB
pub mod build_db;
// pub mod compact_hash;
// pub mod estimate_capacity;
// pub mod kv_store;
// pub mod mmap_file;
// pub mod mmscanner;
// pub mod omp_hack;
// pub mod classify3;
// pub mod aa_translate;
// pub mod seqreader;
// pub mod taxonomy;
// pub mod utilities;

//  lookup accession numbers
// pub mod lookup_accession_numbers;
// pub mod readcounts;
// // pub mod mmap_file;
// // pub mod omp_hack;
// // pub mod utilities;

// k2mask
pub mod k2mask;
pub mod seqreader;
pub mod classify3;
mod compact_hash;
mod estimate_capacity;
mod gz_stream;
mod hyperloglogplus;

fn main() {
    println!("Hello, world!");
    k2mask::main();
    seqreader::main();
}
