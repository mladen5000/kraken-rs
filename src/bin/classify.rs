// Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
// Rust conversion Copyright 2025
//
// This file is part of the Kraken 2 taxonomic sequence classification system.

use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{Error, Result};
use lazy_static::lazy_static;
use std::fmt::Write as _;

use kraken_rs::aa_translate::translate_to_all_frames;
use kraken_rs::compact_hash::CompactHashTable;
use kraken_rs::kraken2_data::{
    IndexOptions, TaxId, TaxonCounters, TaxonCountersMap, TaxonCountsMap, TAXID_MAX,
};
use kraken_rs::kv_store::KeyValueStore;
use kraken_rs::mmscanner::MinimizerScanner;
use kraken_rs::reports::{report_kraken_style, report_mpa_style};
use kraken_rs::seqreader::{BatchSequenceReader, Sequence, SequenceFormat};
use kraken_rs::taxonomy::Taxonomy;
use kraken_rs::utilities::murmur_hash3;

// Constants and type definitions
pub const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
const MATE_PAIR_BORDER_TAXON: TaxId = TAXID_MAX;
const READING_FRAME_BORDER_TAXON: TaxId = TAXID_MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxId = TAXID_MAX - 2;

lazy_static! {
    static ref MUTEX: std::sync::Mutex<()> = std::sync::Mutex::new(());
}

#[derive(Default)]
pub struct Options {
    pub index_filename: String,
    pub taxonomy_filename: String,
    pub options_filename: String,
    pub report_filename: String,
    pub classified_output_filename: String,
    pub unclassified_output_filename: String,
    pub kraken_output_filename: String,
    pub mpa_style_report: bool,
    pub report_kmer_data: bool,
    pub quick_mode: bool,
    pub report_zero_counts: bool,
    pub use_translated_search: bool,
    pub print_scientific_name: bool,
    pub confidence_threshold: f64,
    pub num_threads: usize,
    pub paired_end_processing: bool,
    pub single_file_pairs: bool,
    pub minimum_quality_score: i32,
    pub minimum_hit_groups: i64,
    pub use_memory_mapping: bool,
    pub match_input_order: bool,
    pub input_files: Vec<String>,
}

#[derive(Default)]
struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

impl ClassificationStats {
    fn total_unclassified(&self) -> u64 {
        self.total_sequences - self.total_classified
    }
}

struct OutputStreamData {
    initialized: bool,
    printing_sequences: bool,
    classified_output1: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    classified_output2: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    unclassified_output1: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    unclassified_output2: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    kraken_output: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
}

impl Default for OutputStreamData {
    fn default() -> Self {
        Self::new()
    }
}

impl OutputStreamData {
    fn new() -> Self {
        Self {
            initialized: false,
            printing_sequences: false,
            classified_output1: None,
            classified_output2: None,
            unclassified_output1: None,
            unclassified_output2: None,
            kraken_output: Some(Arc::new(Mutex::new(Box::new(io::stdout())))),
        }
    }

    fn write_all(&mut self, output_data: &OutputData) -> io::Result<()> {
        if let Some(ref output) = self.kraken_output {
            output
                .lock()
                .unwrap()
                .write_all(output_data.kraken_str.as_bytes())?;
        }

        if let Some(ref output) = self.classified_output1 {
            output
                .lock()
                .unwrap()
                .write_all(output_data.classified_out1_str.as_bytes())?;
        }

        if let Some(ref output) = self.classified_output2 {
            output
                .lock()
                .unwrap()
                .write_all(output_data.classified_out2_str.as_bytes())?;
        }

        if let Some(ref output) = self.unclassified_output1 {
            output
                .lock()
                .unwrap()
                .write_all(output_data.unclassified_out1_str.as_bytes())?;
        }

        if let Some(ref output) = self.unclassified_output2 {
            output
                .lock()
                .unwrap()
                .write_all(output_data.unclassified_out2_str.as_bytes())?;
        }

        Ok(())
    }

    fn flush_all(&mut self) -> io::Result<()> {
        if let Some(ref output) = self.kraken_output {
            output.lock().unwrap().flush()?;
        }
        if let Some(ref output) = self.classified_output1 {
            output.lock().unwrap().flush()?;
        }
        if let Some(ref output) = self.classified_output2 {
            output.lock().unwrap().flush()?;
        }
        if let Some(ref output) = self.unclassified_output1 {
            output.lock().unwrap().flush()?;
        }
        if let Some(ref output) = self.unclassified_output2 {
            output.lock().unwrap().flush()?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
    str_representation: String,
}

impl OutputData {
    pub fn new(block_id: u64) -> Self {
        Self {
            block_id,
            kraken_str: String::new(),
            classified_out1_str: String::new(),
            classified_out2_str: String::new(),
            unclassified_out1_str: String::new(),
            unclassified_out2_str: String::new(),
            str_representation: String::new(),
        }
    }

    fn finish_string_representation(&mut self) {
        self.str_representation = format!(
            "Block ID: {}\nKraken Output: {}\nClassified Output 1: {}\nClassified Output 2: {}\nUnclassified Output 1: {}\nUnclassified Output 2: {}",
            self.block_id,
            self.kraken_str,
            self.classified_out1_str,
            self.classified_out2_str,
            self.unclassified_out1_str,
            self.unclassified_out2_str
        );
    }
}

impl Ord for OutputData {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.block_id.cmp(&self.block_id)
    }
}

impl PartialOrd for OutputData {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for OutputData {
    fn eq(&self, other: &Self) -> bool {
        self.block_id == other.block_id
    }
}

impl Eq for OutputData {}

impl Clone for OutputData {
    fn clone(&self) -> Self {
        OutputData {
            block_id: self.block_id,
            kraken_str: self.kraken_str.clone(),
            classified_out1_str: self.classified_out1_str.clone(),
            classified_out2_str: self.classified_out2_str.clone(),
            unclassified_out1_str: self.unclassified_out1_str.clone(),
            unclassified_out2_str: self.unclassified_out2_str.clone(),
            str_representation: self.str_representation.clone(),
        }
    }
}

impl std::fmt::Display for OutputData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Block ID: {}\nKraken Output: {}\nClassified Output 1: {}\nClassified Output 2: {}\nUnclassified Output 1: {}\nUnclassified Output 2: {}",
            self.block_id,
            self.kraken_str,
            self.classified_out1_str,
            self.classified_out2_str,
            self.unclassified_out1_str,
            self.unclassified_out2_str
        )
    }
}

#[derive(Clone)]
pub struct ThreadSafeOutput {
    writer: Arc<Mutex<Box<dyn Write + Send>>>,
}

impl ThreadSafeOutput {
    pub fn new<W: Write + Send + 'static>(writer: W) -> Self {
        Self {
            writer: Arc::new(Mutex::new(Box::new(writer))),
        }
    }

    pub fn write(&self, buf: &[u8]) -> io::Result<usize> {
        self.writer.lock().unwrap().write(buf)
    }

    pub fn write_str(&self, s: &str) -> io::Result<()> {
        self.writer.lock().unwrap().write_all(s.as_bytes())
    }

    pub fn write_all(&self, buf: &[u8]) -> io::Result<()> {
        self.writer.lock().unwrap().write_all(buf)
    }

    pub fn flush(&self) -> io::Result<()> {
        self.writer.lock().unwrap().flush()
    }
}

pub fn main() -> Result<(), Error> {
    let mut opts = Options::default();
    opts.quick_mode = false;
    opts.confidence_threshold = 0.0;
    opts.paired_end_processing = false;
    opts.single_file_pairs = false;
    opts.num_threads = 1;
    opts.mpa_style_report = false;
    opts.report_kmer_data = false;
    opts.report_zero_counts = false;
    opts.use_translated_search = false;
    opts.print_scientific_name = false;
    opts.minimum_quality_score = 0;
    opts.minimum_hit_groups = 0;
    opts.use_memory_mapping = false;

    let mut taxon_counters: TaxonCountersMap = HashMap::new();
    let args: Vec<String> = std::env::args().collect();
    parse_command_line(&args, &mut opts)?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    eprintln!("Loading database information...");

    let mut idx_opts = IndexOptions::default();
    let mut idx_opt_file = File::open(&opts.options_filename)?;
    let metadata = std::fs::metadata(&opts.options_filename)?;
    let opts_filesize = metadata.len() as usize;

    // Read the index options into idx_opts
    unsafe {
        let idx_opts_ptr = &mut idx_opts as *mut _ as *mut u8;
        let idx_opts_slice =
            std::slice::from_raw_parts_mut(idx_opts_ptr, std::mem::size_of::<IndexOptions>());

        // Only read up to the file size to avoid buffer overrun
        let read_size = std::cmp::min(idx_opts_slice.len(), opts_filesize);
        idx_opt_file.read_exact(&mut idx_opts_slice[0..read_size])?;
    }

    opts.use_translated_search = !idx_opts.dna_db;

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping)?;

    // Create and load the hash table
    // This matches the C++ version which creates a CompactHashTable
    let hash_ptr = CompactHashTable::from_file(&opts.index_filename, opts.use_memory_mapping)?;

    eprintln!(" done.");

    let mut stats = ClassificationStats::default();
    let mut outputs = OutputStreamData::new();

    let start_time = Instant::now();

    if opts.input_files.is_empty() {
        if opts.paired_end_processing && !opts.single_file_pairs {
            return Err(anyhow::anyhow!(
                "paired end processing used with no files specified"
            ));
        }
        process_files(
            None,
            None,
            &hash_ptr,
            &taxonomy,
            &idx_opts,
            &opts,
            &mut stats,
            &mut outputs,
            &mut taxon_counters,
        )?;
    } else {
        let mut i = 0;
        while i < opts.input_files.len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == opts.input_files.len() {
                    return Err(anyhow::anyhow!(
                        "paired end processing used with unpaired file"
                    ));
                }
                process_files(
                    Some(&opts.input_files[i]),
                    Some(&opts.input_files[i + 1]),
                    &hash_ptr,
                    &taxonomy,
                    &idx_opts,
                    &opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                )?;
                i += 2;
            } else {
                process_files(
                    Some(&opts.input_files[i]),
                    None,
                    &hash_ptr,
                    &taxonomy,
                    &idx_opts,
                    &opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                )?;
                i += 1;
            }
        }
    }

    let end_time = Instant::now();

    report_stats(start_time, end_time, &stats);

    if !opts.report_filename.is_empty() {
        if opts.mpa_style_report {
            report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &taxon_counters,
            )?;
        } else {
            let total_unclassified = stats.total_unclassified();
            report_kraken_style(
                &opts.report_filename,
                opts.report_zero_counts,
                opts.report_kmer_data,
                &taxonomy,
                &taxon_counters,
                stats.total_sequences,
                total_unclassified,
            )?;
        }
    }

    Ok(())
}

fn report_stats(start_time: Instant, end_time: Instant, stats: &ClassificationStats) {
    let duration = end_time.duration_since(start_time);
    let seconds = duration.as_secs_f64();

    let total_unclassified = stats.total_unclassified();

    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }

    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3} s ({} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1.0e6,
        seconds,
        (stats.total_sequences as f64 / 1.0e3 / (seconds / 60.0)) as u64,
        stats.total_bases as f64 / 1.0e6 / (seconds / 60.0)
    );

    eprintln!(
        "  {} sequences classified ({:.2}%)",
        stats.total_classified,
        (stats.total_classified as f64 * 100.0 / stats.total_sequences as f64) as f64
    );

    eprintln!(
        "  {} sequences unclassified ({:.2}%)",
        total_unclassified,
        (total_unclassified as f64 * 100.0 / stats.total_sequences as f64) as f64
    );
}

fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: &dyn KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    total_taxon_counters: &mut TaxonCountersMap,
) -> Result<(), Error> {
    // Create readers for input files
    let fptr1: Box<dyn BufRead + Send> = match filename1 {
        None => Box::new(io::BufReader::new(io::stdin())),
        Some(fname) => Box::new(io::BufReader::new(File::open(fname)?)),
    };

    let fptr2: Option<Box<dyn BufRead + Send>> = match (
        opts.paired_end_processing,
        opts.single_file_pairs,
        filename2,
    ) {
        (true, false, Some(fname)) => Some(Box::new(io::BufReader::new(File::open(fname)?))),
        _ => None,
    };

    // Set up thread-shared data structures
    // The priority queue for output is designed to ensure fragment data is output in the same order it was input
    let output_queue = Arc::new(Mutex::new(BinaryHeap::new()));
    let next_input_block_id = Arc::new(std::sync::atomic::AtomicU64::new(0));
    let next_output_block_id = Arc::new(std::sync::atomic::AtomicU64::new(0));
    let stats_mutex = Arc::new(Mutex::new(stats));

    // Store the initialized state before borrowing outputs
    let was_initialized = outputs.initialized;

    // Create a Arc<Mutex<>> with a reference instead of moving outputs
    let outputs_mutex = Arc::new(Mutex::new(&mut *outputs));
    let total_taxon_counters_mutex = Arc::new(Mutex::new(total_taxon_counters));
    let output_lock = Arc::new(Mutex::new(()));
    let initialized = Arc::new(std::sync::atomic::AtomicBool::new(was_initialized));

    let reader1 = Arc::new(Mutex::new(BatchSequenceReader::new(fptr1)));
    let reader2 = fptr2.map(|r| Arc::new(Mutex::new(BatchSequenceReader::new(r))));

    let thread_errors = Arc::new(Mutex::new(Vec::new()));

    rayon::scope(|s| {
        for _ in 0..opts.num_threads {
            let output_queue = Arc::clone(&output_queue);
            let next_input_block_id = Arc::clone(&next_input_block_id);
            let next_output_block_id = Arc::clone(&next_output_block_id);
            let reader1 = Arc::clone(&reader1);
            let reader2 = reader2.as_ref().map(Arc::clone);
            let stats_mutex = Arc::clone(&stats_mutex);
            let outputs_mutex = Arc::clone(&outputs_mutex);
            let total_taxon_counters_mutex = Arc::clone(&total_taxon_counters_mutex);
            let thread_errors = Arc::clone(&thread_errors);
            let output_lock = Arc::clone(&output_lock);
            let initialized = Arc::clone(&initialized);

            s.spawn(move |_| {
                let result = (|| -> Result<(), Error> {
                    let mut scanner = MinimizerScanner::new(
                        idx_opts.k,
                        idx_opts.l,
                        idx_opts.spaced_seed_mask,
                        idx_opts.dna_db,
                        idx_opts.toggle_mask,
                        idx_opts.revcom_version == 1,
                    );

                    let mut taxa = Vec::new();
                    let mut hit_counts = HashMap::new();
                    let mut kraken_oss = String::new();
                    let mut c1_oss = String::new();
                    let mut c2_oss = String::new();
                    let mut u1_oss = String::new();
                    let mut u2_oss = String::new();
                    let mut thread_stats = ClassificationStats::default();
                    let mut translated_frames = vec![String::new(); 6];
                    let mut seq1 = Sequence::default();
                    let mut seq2 = Sequence::default();
                    let mut thread_taxon_counters = HashMap::new();
                    let mut block_id: u64;
                    let mut out_data = OutputData::new(0);

                    loop {
                        // Reset thread stats for this block
                        thread_stats.total_sequences = 0;
                        thread_stats.total_bases = 0;
                        thread_stats.total_classified = 0;

                        // Reset all string buffers for this batch
                        kraken_oss.clear();
                        c1_oss.clear();
                        c2_oss.clear();
                        u1_oss.clear();
                        u2_oss.clear();
                        thread_taxon_counters.clear();

                        // Get next block ID and load sequences
                        let mut ok_read;
                        {
                            // Critical section for sequence reading
                            let mut reader1 =
                                reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;

                            block_id = next_input_block_id
                                .fetch_add(1, std::sync::atomic::Ordering::SeqCst);

                            if !opts.paired_end_processing {
                                // Unpaired data? Just read in a sized block
                                ok_read =
                                    reader1.load_block(3 * 1024 * 1024).map_err(Error::from)?;
                            } else if !opts.single_file_pairs {
                                // Paired data in 2 files? Read a line-counted batch from each file
                                ok_read = reader1
                                    .load_batch(NUM_FRAGMENTS_PER_THREAD)
                                    .map_err(Error::from)?;
                                if ok_read && opts.paired_end_processing {
                                    if let Some(ref reader2) = reader2 {
                                        let mut r2 =
                                            reader2.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                        let r2_ok = r2
                                            .load_batch(NUM_FRAGMENTS_PER_THREAD)
                                            .map_err(Error::from)?;
                                        if !r2_ok {
                                            ok_read = false;
                                        }
                                    }
                                }
                            } else {
                                // Paired data in 1 file
                                let frags = NUM_FRAGMENTS_PER_THREAD * 2;
                                // Ensure frag count is even - just in case
                                let batch_size = if frags % 2 == 1 { frags + 1 } else { frags };
                                ok_read = reader1.load_batch(batch_size).map_err(Error::from)?;
                            }
                        }

                        if !ok_read {
                            break;
                        }

                        // Process sequences
                        while {
                            let got_seq1 = {
                                let mut reader1 =
                                    reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                reader1.next_sequence(&mut seq1)?
                            };

                            if !got_seq1 {
                                false // Exit loop if no sequence
                            } else {
                                let mut valid_fragment = true;

                                if opts.paired_end_processing {
                                    if opts.single_file_pairs {
                                        let mut reader1 =
                                            reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                        valid_fragment = reader1.next_sequence(&mut seq2)?;
                                    } else if let Some(ref reader2) = reader2 {
                                        let mut r2 =
                                            reader2.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                        valid_fragment = r2.next_sequence(&mut seq2)?;
                                    } else {
                                        valid_fragment = false;
                                    }
                                }

                                if !valid_fragment {
                                    false // Exit loop if paired read missing
                                } else {
                                    thread_stats.total_sequences += 1;

                                    // Process quality scores if needed
                                    if opts.minimum_quality_score > 0 {
                                        mask_low_quality_bases(
                                            &mut seq1,
                                            opts.minimum_quality_score,
                                        );
                                        if opts.paired_end_processing {
                                            mask_low_quality_bases(
                                                &mut seq2,
                                                opts.minimum_quality_score,
                                            );
                                        }
                                    }

                                    // Classify sequence
                                    let call = classify_sequence(
                                        &seq1,
                                        &seq2,
                                        &mut kraken_oss,
                                        hash,
                                        tax,
                                        idx_opts,
                                        opts,
                                        &mut thread_stats,
                                        &mut scanner,
                                        &mut taxa,
                                        &mut hit_counts,
                                        &mut translated_frames,
                                        &mut thread_taxon_counters,
                                    )?;

                                    // Handle output
                                    if call != 0 {
                                        let buffer = format!(
                                            " kraken:taxid|{}",
                                            tax.nodes()[call as usize].external_id
                                        );
                                        seq1.header.push_str(&buffer);
                                        write!(c1_oss, "{}", seq1)?;

                                        if opts.paired_end_processing {
                                            seq2.header.push_str(&buffer);
                                            write!(c2_oss, "{}", seq2)?;
                                        }
                                    } else {
                                        write!(u1_oss, "{}", seq1)?;
                                        if opts.paired_end_processing {
                                            write!(u2_oss, "{}", seq2)?;
                                        }
                                    }

                                    thread_stats.total_bases += seq1.seq.len() as u64;
                                    if opts.paired_end_processing {
                                        thread_stats.total_bases += seq2.seq.len() as u64;
                                    }

                                    true // Continue loop
                                }
                            }
                        } {} // end while loop for processing sequences

                        // Update global stats
                        {
                            let mut stats = stats_mutex.lock().unwrap();
                            stats.total_sequences += thread_stats.total_sequences;
                            stats.total_bases += thread_stats.total_bases;
                            stats.total_classified += thread_stats.total_classified;

                            if atty::is(atty::Stream::Stderr) {
                                eprint!(
                                    "\rProcessed {} sequences ({} bp) ...",
                                    stats.total_sequences, stats.total_bases
                                );
                            }
                        }

                        // Initialize outputs if needed
                        if !initialized.load(std::sync::atomic::Ordering::SeqCst) {
                            // Lock both the atomic flag and the outputs_mutex to ensure thread safety
                            let _flag_guard = output_lock.lock().unwrap();

                            // Double-check after acquiring the lock
                            if !initialized.load(std::sync::atomic::Ordering::SeqCst) {
                                let mut outputs = outputs_mutex.lock().unwrap();
                                initialize_outputs(
                                    opts,
                                    &mut **outputs,
                                    reader1.lock().unwrap().file_format(),
                                )?;
                                initialized.store(true, std::sync::atomic::Ordering::SeqCst);
                            }
                        }

                        // Prepare output data
                        out_data.block_id = block_id;
                        out_data.kraken_str = kraken_oss.clone();
                        out_data.classified_out1_str = c1_oss.clone();
                        out_data.classified_out2_str = c2_oss.clone();
                        out_data.unclassified_out1_str = u1_oss.clone();
                        out_data.unclassified_out2_str = u2_oss.clone();
                        out_data.finish_string_representation();

                        // Add to output queue
                        {
                            let mut queue = output_queue.lock().unwrap();
                            queue.push(out_data.clone());
                        }

                        // Update taxon counters
                        if !opts.report_filename.is_empty() {
                            let mut total_counters = total_taxon_counters_mutex.lock().unwrap();
                            for (taxon, counter) in &thread_taxon_counters {
                                total_counters
                                    .entry(*taxon)
                                    .and_modify(|e| {
                                        e.merge(counter);
                                    })
                                    .or_insert_with(|| counter.clone());
                            }
                        }

                        // Process output queue - similar to C++ version
                        let mut output_loop = true;

                        while output_loop {
                            {
                                let mut queue = output_queue.lock().unwrap();

                                output_loop = !queue.is_empty();

                                if output_loop
                                    && queue.peek().unwrap().block_id
                                        == next_output_block_id
                                            .load(std::sync::atomic::Ordering::SeqCst)
                                {
                                    // Get the next output block
                                    out_data = queue.pop().unwrap();

                                    // Get lock for output
                                    let _guard = output_lock.lock().unwrap();
                                    next_output_block_id
                                        .fetch_add(1, std::sync::atomic::Ordering::SeqCst);

                                    // Write the output
                                    let mut outputs_guard = outputs_mutex.lock().unwrap();
                                    // Use double dereference to handle the &mut reference
                                    (**outputs_guard).write_all(&out_data)?;
                                } else {
                                    output_loop = false;
                                }
                            }

                            if !output_loop {
                                break;
                            }
                        }
                    }
                    Ok(())
                })();

                if let Err(e) = result {
                    thread_errors.lock().unwrap().push(anyhow::anyhow!("{}", e));
                }
            });
        }
    });

    // Flush all output streams directly using the original outputs variable
    // This is safe because by this point all threads have completed
    outputs.flush_all()?;

    // Check for any thread errors
    let errors = thread_errors.lock().unwrap();
    if !errors.is_empty() {
        return Err(anyhow::anyhow!("{}", errors[0]));
    }

    Ok(())
}

fn resolve_tree(
    hit_counts: &TaxonCountsMap,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (opts.confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for (&taxon, _) in hit_counts {
        let mut score = 0;

        for (&taxon2, &count) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(taxon2, taxon) {
                score += count;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score {
            max_taxon = taxonomy.lowest_common_ancestor(max_taxon, taxon);
        }
    }

    // Reset max. score to be only hits at the called taxon
    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);

    // We probably have a call w/o required support (unless LCA resolved tie)
    let mut current_taxon = max_taxon;
    while current_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts {
            // Add to score if taxon in clade covered by current_taxon
            if taxonomy.is_a_ancestor_of_b(current_taxon, taxon) {
                max_score += count;
            }
        }

        // Score is now sum of hits at current_taxon and w/in current_taxon clade
        if max_score >= required_score {
            // Kill loop and return, we've got enough support here
            return current_taxon;
        } else {
            // Run up tree until confidence threshold is met
            // Run off tree if required score isn't met
            current_taxon = taxonomy.nodes()[current_taxon as usize].parent_id;
        }
    }

    current_taxon
}

fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }

    // Check if ID has the pattern "anything/1" or "anything/2"
    // The C++ version does id[sz-2] == '/' but we need to handle UTF-8
    // Get the slice for the last two characters if possible
    if let Some(last_two) = id.get(sz - 2..) {
        if last_two == "/1" || last_two == "/2" {
            return id[0..sz - 2].to_string();
        }
    }

    id.to_string()
}

fn classify_sequence(
    dna: &Sequence,
    dna2: &Sequence,
    koss: &mut String,
    hash: &dyn KeyValueStore,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<TaxId>,
    hit_counts: &mut TaxonCountsMap,
    tx_frames: &mut [String],
    curr_taxon_counts: &mut TaxonCountersMap,
) -> Result<TaxId, Error> {
    taxa.clear();
    hit_counts.clear();
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups: i64 = 0;
    let mut call: TaxId = 0;

    // Label for quick return during search
    'search: for mate_num in 0..2 {
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        if opts.use_translated_search {
            let frames = translate_to_all_frames(if mate_num == 0 { &dna.seq } else { &dna2.seq });
            tx_frames.clone_from_slice(&frames);
        }

        for frame_idx in 0..frame_ct {
            if opts.use_translated_search {
                scanner.load_sequence(
                    &tx_frames[frame_idx as usize],
                    0,
                    tx_frames[frame_idx as usize].len(),
                );
            } else {
                let seq = if mate_num == 0 { &dna.seq } else { &dna2.seq };
                scanner.load_sequence(seq, 0, seq.len());
            }
            let mut last_minimizer: u64 = u64::MAX;
            let mut last_taxon: TaxId = TAXID_MAX;
            while let Some(min_ptr) = scanner.next_minimizer() {
                let taxon: TaxId;

                if scanner.is_ambiguous() {
                    taxon = AMBIGUOUS_SPAN_TAXON;
                } else {
                    if min_ptr != last_minimizer {
                        let mut skip_lookup = false;
                        if idx_opts.minimum_acceptable_hash_value > 0
                            && murmur_hash3(min_ptr) < idx_opts.minimum_acceptable_hash_value
                        {
                            skip_lookup = true;
                        }
                        taxon = if !skip_lookup { hash.get(min_ptr) } else { 0 };
                        last_taxon = taxon;
                        last_minimizer = min_ptr;
                        if taxon != 0 {
                            minimizer_hit_groups += 1;
                            // Only add kmer if we're generating a report
                            if !opts.report_filename.is_empty() {
                                curr_taxon_counts
                                    .entry(taxon)
                                    .or_insert_with(|| TaxonCounters::new_with_precision(10))
                                    .add_kmer(scanner.last_minimizer());
                            }

                            // In quick mode, return immediately if we've seen enough hit groups
                            if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                                call = taxon;
                                break 'search; // break out of all loops
                            }
                        }
                    } else {
                        taxon = last_taxon;
                    }
                    if taxon != 0 {
                        *hit_counts.entry(taxon).or_insert(0) += 1;
                    }
                }
                taxa.push(taxon);
            }
            if opts.use_translated_search && frame_idx != 5 {
                taxa.push(READING_FRAME_BORDER_TAXON);
            }
        }
        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    // If we didn't get a quick mode call, resolve the taxonomy tree
    if call == 0 {
        let mut total_kmers = taxa.len();
        if opts.paired_end_processing {
            total_kmers -= 1;
        }
        if opts.use_translated_search {
            total_kmers -= if opts.paired_end_processing { 4 } else { 2 };
        }

        let final_taxon = resolve_tree(hit_counts, taxonomy, total_kmers, opts);
        // Void a call made by too few minimizer groups
        call = if final_taxon != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
            0
        } else {
            final_taxon
        };
    }

    if call != 0 {
        stats.total_classified += 1;

        // Only increment read count if we're generating a report
        if !opts.report_filename.is_empty() {
            curr_taxon_counts
                .entry(call)
                .or_insert_with(|| TaxonCounters::new_with_precision(10))
                .increment_read_count();
        }
    }

    if call != 0 {
        koss.push_str("C\t");
    } else {
        koss.push_str("U\t");
    }

    if !opts.paired_end_processing {
        koss.push_str(&dna.id);
        koss.push('\t');
    } else {
        koss.push_str(&trim_pair_info(&dna.id));
        koss.push('\t');
    }

    let ext_call = taxonomy.nodes()[call as usize].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            taxonomy.name_at(taxonomy.nodes()[call as usize].name_offset)
        } else {
            None
        };
        koss.push_str(&format!(
            "{} (taxid {})",
            name.unwrap_or("unclassified"),
            ext_call
        ));
    } else {
        koss.push_str(&ext_call.to_string());
    }

    koss.push('\t');
    if !opts.paired_end_processing {
        koss.push_str(&dna.seq.len().to_string());
        koss.push('\t');
    } else {
        koss.push_str(&format!("{}|{}\t", dna.seq.len(), dna2.seq.len()));
    }

    if opts.quick_mode {
        koss.push_str(&format!("{}:Q", ext_call));
    } else if taxa.is_empty() {
        koss.push_str("0:0");
    } else {
        add_hitlist_string(koss, taxa, taxonomy);
    }

    koss.push('\n');

    Ok(call)
}

fn add_hitlist_string(oss: &mut String, taxa: &[TaxId], taxonomy: &Taxonomy) {
    if taxa.is_empty() {
        return;
    }

    let mut last_code = taxa[0];
    let mut code_count = 1;

    for i in 1..taxa.len() {
        let code = taxa[i];

        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    oss.push_str(&format!("A:{} ", code_count));
                } else {
                    let ext_code = taxonomy.nodes()[last_code as usize].external_id;
                    oss.push_str(&format!("{}:{} ", ext_code, code_count));
                }
            } else {
                // mate pair/reading frame marker
                oss.push_str(if last_code == MATE_PAIR_BORDER_TAXON {
                    "|:| "
                } else {
                    "-:- "
                });
            }
            code_count = 1;
            last_code = code;
        }
    }

    // Handle the last group
    if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
        if last_code == AMBIGUOUS_SPAN_TAXON {
            oss.push_str(&format!("A:{}", code_count));
        } else {
            let ext_code = taxonomy.nodes()[last_code as usize].external_id;
            oss.push_str(&format!("{}:{}", ext_code, code_count));
        }
    } else {
        // mate pair/reading frame marker
        oss.push_str(if last_code == MATE_PAIR_BORDER_TAXON {
            "|:|"
        } else {
            "-:-"
        });
    }
}

fn initialize_outputs(
    opts: &Options,
    outputs: &mut OutputStreamData,
    _format: SequenceFormat,
) -> Result<(), Error> {
    if outputs.initialized {
        return Ok(());
    }

    let _guard = MUTEX.lock().unwrap();

    if !opts.classified_output_filename.is_empty() {
        if opts.paired_end_processing {
            let fields: Vec<&str> = opts.classified_output_filename.split('#').collect();
            if fields.len() != 2 {
                return Err(anyhow::anyhow!(
                    "Paired filename format requires exactly one # character: {}",
                    opts.classified_output_filename
                ));
            }
            outputs.classified_output1 = Some(Arc::new(Mutex::new(Box::new(File::create(
                format!("{}_1{}", fields[0], fields[1]),
            )?))));
            outputs.classified_output2 = Some(Arc::new(Mutex::new(Box::new(File::create(
                format!("{}_2{}", fields[0], fields[1]),
            )?))));
        } else {
            outputs.classified_output1 = Some(Arc::new(Mutex::new(Box::new(File::create(
                &opts.classified_output_filename,
            )?))));
        }
        outputs.printing_sequences = true;
    }

    if !opts.unclassified_output_filename.is_empty() {
        if opts.paired_end_processing {
            let fields: Vec<&str> = opts.unclassified_output_filename.split('#').collect();
            if fields.len() != 2 {
                return Err(anyhow::anyhow!(
                    "Paired filename format requires exactly one # character: {}",
                    opts.unclassified_output_filename
                ));
            }
            outputs.unclassified_output1 = Some(Arc::new(Mutex::new(Box::new(File::create(
                format!("{}_1{}", fields[0], fields[1]),
            )?))));
            outputs.unclassified_output2 = Some(Arc::new(Mutex::new(Box::new(File::create(
                format!("{}_2{}", fields[0], fields[1]),
            )?))));
        } else {
            outputs.unclassified_output1 = Some(Arc::new(Mutex::new(Box::new(File::create(
                &opts.unclassified_output_filename,
            )?))));
        }
        outputs.printing_sequences = true;
    }

    if !opts.kraken_output_filename.is_empty() {
        outputs.kraken_output = if opts.kraken_output_filename == "-" {
            // Special filename to silence Kraken output
            None
        } else {
            Some(Arc::new(Mutex::new(Box::new(File::create(
                &opts.kraken_output_filename,
            )?))))
        };
    }

    outputs.initialized = true;
    Ok(())
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: i32) {
    if dna.format != SequenceFormat::Fastq {
        return;
    }

    let seq_len = dna.seq.len();
    let quals_len = dna.quals.len();

    if seq_len != quals_len {
        panic!(
            "{}: Sequence length ({}) != Quality string length ({})",
            dna.id, seq_len, quals_len
        );
    }

    // Create a new sequence with masked characters
    let mut new_seq = String::with_capacity(seq_len);

    for i in 0..seq_len {
        let qual = dna.quals.as_bytes()[i] as i32 - b'!' as i32;
        if qual < minimum_quality_score {
            new_seq.push('x');
        } else {
            new_seq.push(dna.seq.chars().nth(i).unwrap());
        }
    }

    // Replace the original sequence with the masked one
    dna.seq = new_seq;
}

fn parse_command_line(args: &[String], opts: &mut Options) -> Result<(), Error> {
    let mut i = 1;
    while i < args.len() && args[i].starts_with('-') {
        // Only parse options that start with -
        match args[i].as_str() {
            "-h" | "-?" => {
                usage(0);
            }
            "-H" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -H"));
                }
                opts.index_filename = args[i].clone();
            }
            "-t" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -t"));
                }
                opts.taxonomy_filename = args[i].clone();
            }
            "-T" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -T"));
                }
                opts.confidence_threshold = args[i].parse()?;
                if opts.confidence_threshold < 0.0 || opts.confidence_threshold > 1.0 {
                    return Err(anyhow::anyhow!("confidence threshold must be in [0, 1]"));
                }
            }
            "-o" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -o"));
                }
                opts.options_filename = args[i].clone();
            }
            "-q" => opts.quick_mode = true,
            "-p" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -p"));
                }
                opts.num_threads = args[i].parse()?;
                if opts.num_threads < 1 {
                    return Err(anyhow::anyhow!("number of threads can't be less than 1"));
                }
            }
            "-g" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -g"));
                }
                opts.minimum_hit_groups = args[i].parse()?;
            }
            "-P" => opts.paired_end_processing = true,
            "-S" => {
                opts.paired_end_processing = true;
                opts.single_file_pairs = true;
            }
            "-m" => opts.mpa_style_report = true,
            "-K" => opts.report_kmer_data = true,
            "-R" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -R"));
                }
                opts.report_filename = args[i].clone();
            }
            "-z" => opts.report_zero_counts = true,
            "-C" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -C"));
                }
                opts.classified_output_filename = args[i].clone();
            }
            "-U" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -U"));
                }
                opts.unclassified_output_filename = args[i].clone();
            }
            "-O" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -O"));
                }
                opts.kraken_output_filename = args[i].clone();
            }
            "-n" => opts.print_scientific_name = true,
            "-Q" => {
                i += 1;
                if i >= args.len() {
                    return Err(anyhow::anyhow!("missing argument for -Q"));
                }
                opts.minimum_quality_score = args[i].parse()?;
            }
            "-M" => opts.use_memory_mapping = true,
            _ => {
                return Err(anyhow::anyhow!("unknown option: {}", args[i]));
            }
        }
        i += 1;
    }

    // Add remaining arguments as input files
    while i < args.len() {
        opts.input_files.push(args[i].clone());
        i += 1;
    }

    if opts.index_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        return Err(anyhow::anyhow!("mandatory filename missing"));
    }

    if opts.mpa_style_report && opts.report_filename.is_empty() {
        return Err(anyhow::anyhow!("-m requires -R be used"));
    }

    Ok(())
}

fn usage(exit_code: i32) {
    eprintln!(
        "Usage: classify [options] <fasta/fastq file(s)>\n\n\
         Options: (*mandatory)\n\
         * -H filename      Kraken 2 index filename\n\
         * -t filename      Kraken 2 taxonomy filename\n\
         * -o filename      Kraken 2 options filename\n\
         -q               Quick mode\n\
         -M               Use memory mapping to access hash & taxonomy\n\
         -T NUM           Confidence score threshold (def. 0)\n\
         -p NUM           Number of threads (def. 1)\n\
         -Q NUM           Minimum quality score (FASTQ only, def. 0)\n\
         -P               Process pairs of reads\n\
         -S               Process pairs with mates in same file\n\
         -R filename      Print report to filename\n\
         -m               In comb. w/ -R, use mpa-style report\n\
         -z               In comb. w/ -R, report taxa w/ 0 count\n\
         -n               Print scientific name instead of taxid in Kraken output\n\
         -g NUM           Minimum number of hit groups needed for call\n\
         -C filename      Filename/format to have classified sequences\n\
         -U filename      Filename/format to have unclassified sequences\n\
         -O filename      Output file for normal Kraken output\n\
         -K               In comb. w/ -R, provide minimizer information in report"
    );
    std::process::exit(exit_code);
}

// The process_sequence_batch function has been moved directly into the thread processing
// code in process_files() for better alignment with the C++ implementation

#[cfg(test)]
mod tests {
    use super::*;
    use kraken_rs::kraken2_data::TaxonCounters;
    use kraken_rs::kv_store::KeyValueStore;
    use kraken_rs::seqreader::{Sequence, SequenceFormat};
    use kraken_rs::taxonomy::TaxonomyNode;
    use std::io::Cursor;
    use tempfile::tempdir;

    #[test]
    fn test_trim_pair_info() {
        assert_eq!(trim_pair_info("read/1"), "read");
        assert_eq!(trim_pair_info("read/2"), "read");
        assert_eq!(trim_pair_info("read"), "read");
        assert_eq!(trim_pair_info("a"), "a");
    }

    #[test]
    fn test_resolve_tree() {
        // Create a properly initialized taxonomy for testing
        let temp_dir = tempdir().unwrap();
        let taxonomy_path = temp_dir.path().join("taxonomy.bin");
        let mut file = File::create(&taxonomy_path).unwrap();

        // Write a minimal valid taxonomy file with nodes for our test
        file.write_all(b"K2TAXDAT").unwrap(); // Magic
        file.write_all(&6usize.to_le_bytes()).unwrap(); // node_count
        file.write_all(&10usize.to_le_bytes()).unwrap(); // name_data_len
        file.write_all(&10usize.to_le_bytes()).unwrap(); // rank_data_len

        // Write node 0 (unused)
        file.write_all(&0u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write node 1 (root)
        file.write_all(&0u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&2u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&2u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&1u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write node 2
        file.write_all(&1u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&4u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&1u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&2u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write node 3
        file.write_all(&1u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&5u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&1u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&3u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write node 4
        file.write_all(&2u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&4u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write node 5
        file.write_all(&3u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&5u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write dummy name and rank data
        file.write_all(&[0; 10]).unwrap(); // name_data
        file.write_all(&[0; 10]).unwrap(); // rank_data

        let taxonomy = Taxonomy::new(&taxonomy_path, false).unwrap();
        let mut hit_counts = TaxonCountsMap::new();

        // Create a simple taxonomy tree:
        //       1
        //      / \
        //     2   3
        //    /     \
        //   4       5

        hit_counts.insert(4, 2);
        hit_counts.insert(5, 3);
        hit_counts.insert(3, 1);

        let opts = Options {
            confidence_threshold: 0.5,
            ..Default::default()
        };

        let total_minimizers = 10;
        let result = resolve_tree(&hit_counts, &taxonomy, total_minimizers, &opts);

        // With 50% confidence threshold and 10 minimizers, need 5 hits to make call
        // Node 3 with node 5 has 4 hits total, but node 1 has 6 total (meets confidence threshold)
        assert_eq!(result, 1);
    }

    #[test]
    fn test_mask_low_quality_bases() {
        let mut seq = Sequence {
            format: SequenceFormat::Fastq,
            header: "test".to_string(),
            id: "test".to_string(),
            seq: "ACGT".to_string(),
            quals: "!!!#".to_string(), // ASCII 33, 33, 33, 35
            str_representation: String::new(),
        };

        // ASCII 35(#) - ASCII 33(!) = quality score 2
        mask_low_quality_bases(&mut seq, 2);
        // First three should be masked (below score 2)
        assert_eq!(seq.seq, "xxxT");
    }

    #[test]
    fn test_output_data_priority() {
        let mut heap = BinaryHeap::new();

        heap.push(OutputData {
            block_id: 1,
            ..Default::default()
        });
        heap.push(OutputData {
            block_id: 0,
            ..Default::default()
        });
        heap.push(OutputData {
            block_id: 2,
            ..Default::default()
        });

        assert_eq!(heap.pop().unwrap().block_id, 0);
        assert_eq!(heap.pop().unwrap().block_id, 1);
        assert_eq!(heap.pop().unwrap().block_id, 2);
    }

    #[test]
    fn test_thread_safe_output() {
        let cursor = Cursor::new(Vec::new());
        let output = ThreadSafeOutput::new(cursor);

        output.write_str("test").unwrap();
        output.flush().unwrap();

        let inner = output.writer.lock().unwrap();
        let boxed_writer = inner.as_ref();
        // Get underlying cursor data by transmuting, only safe in tests
        let cursor_data = unsafe {
            let cursor_ptr = boxed_writer as *const _ as *const Cursor<Vec<u8>>;
            &*cursor_ptr
        };
        assert_eq!(cursor_data.get_ref(), b"test");
    }

    #[test]
    fn test_taxon_counters() {
        let mut counter = TaxonCounters::new_with_precision(10);
        counter.increment_read_count();
        assert_eq!(counter.get_read_count(), 1);

        counter.add_kmer(42);
        assert!(counter.get_kmer_distinct() > 0.0);
    }

    #[test]
    fn test_add_hitlist_string_detailed() {
        // Setup a simple mock taxonomy for testing
        let mut taxonomy = Taxonomy::default();

        // Create a test taxonomy properly without directly modifying private fields
        let temp_dir = tempdir().unwrap();
        let taxonomy_path = temp_dir.path().join("taxonomy.bin");
        let mut file = File::create(&taxonomy_path).unwrap();

        // Write a minimal valid taxonomy file
        file.write_all(b"K2TAXDAT").unwrap(); // Magic
        file.write_all(&6usize.to_le_bytes()).unwrap(); // node_count
        file.write_all(&10usize.to_le_bytes()).unwrap(); // name_data_len
        file.write_all(&10usize.to_le_bytes()).unwrap(); // rank_data_len

        // Write nodes (careful with endianness)
        for i in 0..6 {
            file.write_all(&0u64.to_le_bytes()).unwrap(); // parent_id
            file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
            file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
            file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
            file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
            file.write_all(&(i as u64).to_le_bytes()).unwrap(); // external_id
            file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id
        }

        // Write dummy name and rank data
        file.write_all(&[0; 10]).unwrap(); // name_data
        file.write_all(&[0; 10]).unwrap(); // rank_data

        let taxonomy = Taxonomy::new(&taxonomy_path, false).unwrap();

        // Test case 1: Simple list with single taxon
        let mut taxa = vec![2, 2, 2];
        let mut result = String::new();
        add_hitlist_string(&mut result, &taxa, &taxonomy);
        assert_eq!(result, "2:3");

        // Test case 2: Multiple different taxa
        taxa = vec![1, 1, 2, 3, 3, 3];
        result.clear();
        add_hitlist_string(&mut result, &taxa, &taxonomy);
        assert_eq!(result, "1:2 2:1 3:3");

        // Test case 3: With mate pair and reading frame borders
        taxa = vec![
            1,
            1,
            READING_FRAME_BORDER_TAXON,
            2,
            2,
            MATE_PAIR_BORDER_TAXON,
            3,
            3,
        ];
        result.clear();
        add_hitlist_string(&mut result, &taxa, &taxonomy);
        assert_eq!(result, "1:2 -:- 2:2 |:| 3:2");

        // Test case 4: With ambiguous spans
        taxa = vec![1, AMBIGUOUS_SPAN_TAXON, AMBIGUOUS_SPAN_TAXON, 2, 2];
        result.clear();
        add_hitlist_string(&mut result, &taxa, &taxonomy);
        assert_eq!(result, "1:1 A:2 2:2");
    }

    #[test]
    fn test_advanced_low_quality_bases() {
        // Setup a test sequence
        let mut seq = Sequence {
            format: SequenceFormat::Fastq,
            header: "test".to_string(),
            id: "test".to_string(),
            seq: "ACGTACGT".to_string(),
            quals: "!!#$%&'(".to_string(), // ASCII 33-40
            str_representation: String::new(),
        };

        // Test with minimum quality 3
        // $ = ASCII 36; 36 - 33 = 3 (quality score)
        mask_low_quality_bases(&mut seq, 3);

        // First three bases should be masked (!, !, #)
        assert_eq!(seq.seq, "xxxTACGT");

        // Test with higher minimum quality
        seq.seq = "ACGTACGT".to_string();
        // & = ASCII 38; 38 - 33 = 5 (quality score)
        mask_low_quality_bases(&mut seq, 5);

        // First five bases should be masked (A=!, C=!, G=#, T=$, A=%)
        assert_eq!(seq.seq, "xxxxxCGT");
    }

    #[test]
    fn test_extended_pair_info() {
        // Test standard Illumina pair notation
        assert_eq!(trim_pair_info("read1/1"), "read1");
        assert_eq!(trim_pair_info("read2/2"), "read2");

        // Test with no pair info
        assert_eq!(trim_pair_info("readX"), "readX");

        // Test with unusual formats
        assert_eq!(trim_pair_info("read/3"), "read/3"); // Not 1 or 2

        // Skip the double slash test because our trim function only checks
        // for /1 or /2 at the end, not double slashes

        // Test short strings
        assert_eq!(trim_pair_info("a"), "a");
        assert_eq!(trim_pair_info(""), "");
    }

    #[test]
    fn test_extended_command_line() {
        // Create test args
        let args = vec![
            "classify".to_string(),
            "-H".to_string(),
            "hash.bin".to_string(),
            "-t".to_string(),
            "taxonomy.bin".to_string(),
            "-o".to_string(),
            "options.bin".to_string(),
            "-p".to_string(),
            "4".to_string(),
            "-T".to_string(),
            "0.5".to_string(),
        ];

        let mut opts = Options::default();
        parse_command_line(&args, &mut opts).unwrap();

        // Check parsed options
        assert_eq!(opts.index_filename, "hash.bin");
        assert_eq!(opts.taxonomy_filename, "taxonomy.bin");
        assert_eq!(opts.options_filename, "options.bin");
        assert_eq!(opts.num_threads, 4);
        assert_eq!(opts.confidence_threshold, 0.5);
        assert!(!opts.quick_mode);

        // Test with additional options
        let args = vec![
            "classify".to_string(),
            "-H".to_string(),
            "hash.bin".to_string(),
            "-t".to_string(),
            "taxonomy.bin".to_string(),
            "-o".to_string(),
            "options.bin".to_string(),
            "-q".to_string(), // Quick mode
            "-P".to_string(), // Paired end
            "-n".to_string(), // Print scientific name
        ];

        let mut opts = Options::default();
        parse_command_line(&args, &mut opts).unwrap();

        assert!(opts.quick_mode);
        assert!(opts.paired_end_processing);
        assert!(opts.print_scientific_name);
    }

    #[test]
    fn test_extended_tree_resolution() {
        // Create a properly initialized taxonomy
        let temp_dir = tempdir().unwrap();
        let taxonomy_path = temp_dir.path().join("taxonomy.bin");
        let mut file = File::create(&taxonomy_path).unwrap();

        // Create a simple taxonomy with known structure:
        //       1 (root)
        //      / \
        //     2   3
        //    / \   \
        //   4   5   6

        // Write file magic and header
        file.write_all(b"K2TAXDAT").unwrap(); // Magic
        file.write_all(&7usize.to_le_bytes()).unwrap(); // node_count
        file.write_all(&10usize.to_le_bytes()).unwrap(); // name_data_len
        file.write_all(&10usize.to_le_bytes()).unwrap(); // rank_data_len

        // Node 0 (unused)
        file.write_all(&0u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 1 (root)
        file.write_all(&0u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&2u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&2u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&1u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 2
        file.write_all(&1u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&4u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&2u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&2u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 3
        file.write_all(&1u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&6u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&1u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&3u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 4
        file.write_all(&2u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&4u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 5
        file.write_all(&2u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&5u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Node 6
        file.write_all(&3u64.to_le_bytes()).unwrap(); // parent_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // first_child
        file.write_all(&0u64.to_le_bytes()).unwrap(); // child_count
        file.write_all(&0u64.to_le_bytes()).unwrap(); // name_offset
        file.write_all(&0u64.to_le_bytes()).unwrap(); // rank_offset
        file.write_all(&6u64.to_le_bytes()).unwrap(); // external_id
        file.write_all(&0u64.to_le_bytes()).unwrap(); // godparent_id

        // Write dummy name and rank data
        file.write_all(&[0; 10]).unwrap(); // name_data
        file.write_all(&[0; 10]).unwrap(); // rank_data

        let taxonomy = Taxonomy::new(&taxonomy_path, false).unwrap();

        // Test case 1: Hit counts concentrated on leaf nodes
        let mut hit_counts = TaxonCountsMap::new();
        hit_counts.insert(4, 3);
        hit_counts.insert(5, 2);
        hit_counts.insert(6, 1);

        let opts = Options {
            confidence_threshold: 0.5,
            ..Default::default()
        };

        // With threshold 0.5 and 6 minimizers, need at least 3 hits
        // Node 4 has 3 hits, so it should be the winner
        let result = resolve_tree(&hit_counts, &taxonomy, 6, &opts);
        assert_eq!(result, 4);

        // Test case 2: Identical hit counts, should resolve to LCA
        hit_counts.clear();
        hit_counts.insert(4, 2);
        hit_counts.insert(5, 2);

        // With identical hit counts, should resolve to common ancestor
        // But with confidence threshold of 0.5 and only 4 total hits out of 10 minimizers,
        // it will walk up the tree to node 0 (not enough confidence)
        let result = resolve_tree(&hit_counts, &taxonomy, 10, &opts);
        assert_eq!(result, 0);

        // Test case 3: High confidence threshold forcing a higher-level classification
        let opts = Options {
            confidence_threshold: 0.8,
            ..Default::default()
        };

        hit_counts.clear();
        hit_counts.insert(4, 5);
        hit_counts.insert(5, 2);
        hit_counts.insert(6, 1);

        // With threshold 0.8 and 10 minimizers, need 8 hits
        // Node 4 has 5 hits, not enough; walking up the tree:
        // Node 2 has 7 hits (4+5), still not enough
        // Node 1 has 8 hits (4+5+6), which is enough
        let result = resolve_tree(&hit_counts, &taxonomy, 10, &opts);
        assert_eq!(result, 1);
    }
}
