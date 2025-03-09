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

use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{
    IndexOptions, TaxId, TaxonCounters, TaxonCountersMap, TaxonCountsMap, TAXID_MAX,
};
use crate::kv_store::KeyValueStore;
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{BatchSequenceReader, Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy;
use crate::utilities::murmur_hash3;

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
            crate::reports::report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &taxon_counters,
            )?;
        } else {
            let total_unclassified = stats.total_unclassified();
            crate::reports::report_kraken_style(
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

    let output_queue = Arc::new(Mutex::new(BinaryHeap::new()));
    let next_input_block_id = Arc::new(std::sync::atomic::AtomicU64::new(0));
    let next_output_block_id = Arc::new(std::sync::atomic::AtomicU64::new(0));
    let stats_mutex = Arc::new(Mutex::new(stats));
    let outputs_mutex = Arc::new(Mutex::new(outputs));
    let total_taxon_counters_mutex = Arc::new(Mutex::new(total_taxon_counters));

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

            s.spawn(move |_| {
                let result = (|| -> Result<(), Error> {
                    let mut scanner = MinimizerScanner::new(
                        idx_opts.k,
                        idx_opts.l,
                        idx_opts.spaced_seed_mask,
                        idx_opts.dna_db, // Pass the same value from options
                        idx_opts.toggle_mask,
                        idx_opts.revcom_version == 1, // Convert u32 to bool
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

                    loop {
                        // Reset thread stats for this block
                        thread_stats.total_sequences = 0;
                        thread_stats.total_bases = 0;
                        thread_stats.total_classified = 0;

                        let block_id =
                            next_input_block_id.fetch_add(1, std::sync::atomic::Ordering::SeqCst);

                        let ok_read = {
                            let mut reader1 =
                                reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                            if !opts.paired_end_processing {
                                reader1.load_block(3 * 1024 * 1024).map_err(Error::from)
                            } else if !opts.single_file_pairs {
                                match reader1
                                    .load_batch(NUM_FRAGMENTS_PER_THREAD)
                                    .map_err(Error::from)?
                                {
                                    true => {
                                        if let Some(ref reader2) = reader2 {
                                            reader2
                                                .lock()
                                                .map_err(|e| anyhow::anyhow!("{}", e))?
                                                .load_batch(NUM_FRAGMENTS_PER_THREAD)
                                                .map_err(Error::from)
                                        } else {
                                            Ok(true)
                                        }
                                    }
                                    false => Ok(false),
                                }
                            } else {
                                let frags = NUM_FRAGMENTS_PER_THREAD * 2;
                                let batch_size = if frags % 2 == 1 { frags + 1 } else { frags };
                                reader1.load_batch(batch_size).map_err(Error::from)
                            }
                        }?;

                        if !ok_read {
                            break;
                        }

                        // Process sequences
                        loop {
                            let got_seq1 = {
                                let mut reader1 =
                                    reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                reader1.next_sequence(&mut seq1)?
                            };

                            if !got_seq1 {
                                break;
                            }

                            let got_seq2 = if opts.paired_end_processing {
                                if opts.single_file_pairs {
                                    let mut reader1 =
                                        reader1.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                    reader1.next_sequence(&mut seq2)?
                                } else if let Some(ref reader2) = reader2 {
                                    let mut reader2 =
                                        reader2.lock().map_err(|e| anyhow::anyhow!("{}", e))?;
                                    reader2.next_sequence(&mut seq2)?
                                } else {
                                    false
                                }
                            } else {
                                true
                            };

                            if opts.paired_end_processing && !got_seq2 {
                                break;
                            }

                            thread_stats.total_sequences += 1;

                            process_sequence_batch(
                                &mut scanner,
                                &mut taxa,
                                &mut hit_counts,
                                &mut kraken_oss,
                                &mut c1_oss,
                                &mut c2_oss,
                                &mut u1_oss,
                                &mut u2_oss,
                                &mut thread_stats,
                                &mut translated_frames,
                                &mut seq1,
                                &mut seq2,
                                &mut thread_taxon_counters,
                                hash,
                                tax,
                                idx_opts,
                                opts,
                            )?;
                        }

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
                        {
                            let mut outputs = outputs_mutex.lock().unwrap();
                            if !outputs.initialized {
                                initialize_outputs(
                                    opts,
                                    &mut *outputs,
                                    reader1.lock().unwrap().file_format(),
                                )?;
                            }
                        }

                        // Queue output data
                        let mut out_data = OutputData {
                            block_id,
                            kraken_str: kraken_oss.clone(),
                            classified_out1_str: c1_oss.clone(),
                            classified_out2_str: c2_oss.clone(),
                            unclassified_out1_str: u1_oss.clone(),
                            unclassified_out2_str: u2_oss.clone(),
                            str_representation: String::new(),
                        };
                        out_data.finish_string_representation();

                        // Reset all dynamically-growing things for next batch
                        kraken_oss.clear();
                        c1_oss.clear();
                        c2_oss.clear();
                        u1_oss.clear();
                        u2_oss.clear();
                        thread_taxon_counters.clear();

                        output_queue.lock().unwrap().push(out_data);

                        // Update taxon counters
                        if !opts.report_filename.is_empty() {
                            let mut total_counters = total_taxon_counters_mutex.lock().unwrap();
                            for (taxon, counter) in thread_taxon_counters.drain() {
                                total_counters
                                    .entry(taxon)
                                    .and_modify(|e| {
                                        e.merge(&counter);
                                    })
                                    .or_insert_with(|| counter);
                            }
                        }

                        // Process output queue
                        loop {
                            let next_output_id =
                                next_output_block_id.load(std::sync::atomic::Ordering::SeqCst);
                            let mut queue = output_queue.lock().unwrap();

                            if let Some(out_data) = queue.peek() {
                                if out_data.block_id != next_output_id {
                                    break;
                                }

                                let out_data = queue.pop().unwrap();
                                next_output_block_id
                                    .fetch_add(1, std::sync::atomic::Ordering::SeqCst);

                                // Write output
                                let mut outputs = outputs_mutex.lock().unwrap();
                                outputs.write_all(&out_data)?;
                            } else {
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
            // Add to score if taxon in max_taxon's clade
            if taxonomy.is_a_ancestor_of_b(current_taxon, taxon) {
                max_score += count;
            }
        }

        // Score is now sum of hits at max_taxon and w/in max_taxon clade
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

    // Optimize to avoid collecting all chars - just check the relevant positions
    if id.chars().nth(sz - 2) == Some('/') {
        let last_char = id.chars().nth(sz - 1);
        if last_char == Some('1') || last_char == Some('2') {
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

    for mate_num in 0..2 {
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
                        }
                    } else {
                        taxon = last_taxon;
                    }
                    if taxon != 0 {
                        if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                            return Ok(taxon);
                        }
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

    let mut total_kmers = taxa.len();
    if opts.paired_end_processing {
        total_kmers -= 1;
    }
    if opts.use_translated_search {
        total_kmers -= if opts.paired_end_processing { 4 } else { 2 };
    }

    let final_taxon = resolve_tree(hit_counts, taxonomy, total_kmers, opts);
    let result = if final_taxon != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        0
    } else {
        final_taxon
    };

    if result != 0 {
        stats.total_classified += 1;

        // Only increment read count if we're generating a report
        if !opts.report_filename.is_empty() {
            curr_taxon_counts
                .entry(result)
                .or_insert_with(|| TaxonCounters::new_with_precision(10))
                .increment_read_count();
        }
    }

    if result != 0 {
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

    let ext_call = taxonomy.nodes()[result as usize].external_id;
    if opts.print_scientific_name {
        let name = if result != 0 {
            taxonomy.name_at(taxonomy.nodes()[result as usize].name_offset)
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
    } else {
        if taxa.is_empty() {
            koss.push_str("0:0");
        } else {
            add_hitlist_string(koss, taxa, taxonomy);
        }
    }

    koss.push('\n');

    Ok(result)
}

fn add_hitlist_string(oss: &mut String, taxa: &[TaxId], taxonomy: &Taxonomy) {
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
    if dna.seq.len() != dna.quals.len() {
        panic!(
            "{}: Sequence length ({}) != Quality string length ({})",
            dna.id,
            dna.seq.len(),
            dna.quals.len()
        );
    }

    for i in 0..dna.seq.len() {
        if (dna.quals.as_bytes()[i] as i32 - b'!' as i32) < minimum_quality_score {
            dna.seq.replace_range(i..i + 1, "x");
        }
    }
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

/// Process a batch of sequences
///
/// This function processes a batch of sequences, classifying each one and updating
/// the appropriate output buffers and statistics.
fn process_sequence_batch(
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<TaxId>,
    hit_counts: &mut TaxonCountsMap,
    kraken_oss: &mut String,
    c1_oss: &mut String,
    c2_oss: &mut String,
    u1_oss: &mut String,
    u2_oss: &mut String,
    thread_stats: &mut ClassificationStats,
    translated_frames: &mut [String],
    seq1: &mut Sequence,
    seq2: &mut Sequence,
    thread_taxon_counters: &mut TaxonCountersMap,
    hash: &dyn KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
) -> Result<(), Error> {
    // Process quality scores if needed
    if opts.minimum_quality_score > 0 {
        mask_low_quality_bases(seq1, opts.minimum_quality_score);
        if opts.paired_end_processing {
            mask_low_quality_bases(seq2, opts.minimum_quality_score);
        }
    }

    // Classify sequence
    let call = classify_sequence(
        seq1,
        seq2,
        kraken_oss,
        hash,
        tax,
        idx_opts,
        opts,
        thread_stats,
        scanner,
        taxa,
        hit_counts,
        translated_frames,
        thread_taxon_counters,
    )?;

    // Handle output
    if call != 0 {
        let buffer = format!(" kraken:taxid|{}", tax.nodes()[call as usize].external_id);
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

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kraken2_data::TaxonCounters;
    use crate::seqreader::{Sequence, SequenceFormat};
    use std::io::Cursor;

    #[test]
    fn test_trim_pair_info() {
        assert_eq!(trim_pair_info("read/1"), "read");
        assert_eq!(trim_pair_info("read/2"), "read");
        assert_eq!(trim_pair_info("read"), "read");
        assert_eq!(trim_pair_info("a"), "a");
    }

    #[test]
    fn test_resolve_tree() {
        let taxonomy = Taxonomy::default();
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
        // Should resolve to taxon 3 since it has most hits in its clade (4)
        assert_eq!(result, 3);
    }

    #[test]
    fn test_mask_low_quality_bases() {
        let mut seq = Sequence {
            format: SequenceFormat::Fastq,
            header: "test".to_string(),
            id: "test".to_string(),
            seq: "ACGT".to_string(),
            quals: "!!!#".to_string(),
            str_representation: String::new(),
        };

        mask_low_quality_bases(&mut seq, 35); // ASCII 35 = '#'
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
}
