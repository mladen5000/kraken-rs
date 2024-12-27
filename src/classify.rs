// lib.rs or main.rs

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use clap::{App, Arg};
use rayon::prelude::*;

use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{taxid_t, TAXID_MAX};
use crate::kv_store::KeyValueStore;
use crate::mmscanner::MinimizerScanner;
use crate::readcounts::TaxonCount;
use crate::reports::{report_kraken_style, report_mpa_style};
use crate::seqreader::{BatchequenceReader, Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy;
use crate::utilities::{MaskLowQualityBases, MurmurHash3, SplitString};

const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
const MATE_PAIR_BORDER_TAXON: taxid_t = TAXID_MAX;
const READING_FRAME_BORDER_TAXON: taxid_t = TAXID_MAX - 1;
const AMBIGUOUS_SPAN_TAXON: taxid_t = TAXID_MAX - 2;

#[derive(Default)]
struct Options {
    index_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    report_filename: String,
    classified_output_filename: String,
    unclassified_output_filename: String,
    kraken_output_filename: String,
    mpa_style_report: bool,
    report_kmer_data: bool,
    quick_mode: bool,
    report_zero_counts: bool,
    use_translated_search: bool,
    print_scientific_name: bool,
    confidence_threshold: f64,
    num_threads: usize,
    paired_end_processing: bool,
    single_file_pairs: bool,
    minimum_quality_score: i32,
    minimum_hit_groups: i32,
    use_memory_mapping: bool,
    match_input_order: bool, // Not used in original but declared for completeness
}

#[derive(Default)]
struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

struct OutputStreamData {
    initialized: bool,
    printing_sequences: bool,
    classified_output1: Option<BufWriter<File>>,
    classified_output2: Option<BufWriter<File>>,
    unclassified_output1: Option<BufWriter<File>>,
    unclassified_output2: Option<BufWriter<File>>,
    kraken_output: Option<BufWriter<File>>,
}

impl Default for OutputStreamData {
    fn default() -> Self {
        Self {
            initialized: false,
            printing_sequences: false,
            classified_output1: None,
            classified_output2: None,
            unclassified_output1: None,
            unclassified_output2: None,
            kraken_output: None,
        }
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
}

struct IndexOptions {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    dna_db: bool,
    toggle_mask: u64,
    revcom_version: u8,
    minimum_acceptable_hash_value: u64,
}

// Dummy implementations for external functions that were in original code
fn report_mpa_style(
    filename: &str,
    report_zero_counts: bool,
    taxonomy: &Taxonomy,
    taxon_counters: &TaxonCountersT,
) -> Result<(), Box<dyn Error>> {
    reports::report_mpa_style(filename, report_zero_counts, taxonomy, taxon_counters)
}

fn ReportKrakenStyle(
    filename: &str,
    report_zero_counts: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    taxon_counters: &TaxonCountersT,
    total_sequences: u64,
    total_unclassified: u64,
) -> Result<(), Box<dyn Error>> {
    reports::report_kraken_style(
        filename,
        report_zero_counts,
        report_kmer_data,
        taxonomy,
        taxon_counters,
        total_sequences,
        total_unclassified,
    )
}

fn initialize_outputs(
    opts: &Options,
    outputs: &mut OutputStreamData,
    format: SequenceFormat,
) -> Result<(), Box<dyn Error>> {
    if (!outputs.initialized) {
        if (!opts.classified_output_filename.is_empty()) {
            if (opts.paired_end_processing) {
                let fields = SplitString(&opts.classified_output_filename, "#");
                if (fields.len() != 2) {
                    return Err("Paired filename format missing or extra '#'".into());
                }
                outputs.classified_output1 = Some(BufWriter::new(File::create(format!(
                    "{}_1{}",
                    fields[0], fields[1]
                ))?));
                outputs.classified_output2 = Some(BufWriter::new(File::create(format!(
                    "{}_2{}",
                    fields[0], fields[1]
                ))?));
            } else {
                outputs.classified_output1 = Some(BufWriter::new(File::create(
                    &opts.classified_output_filename,
                )?));
            }
            outputs.printing_sequences = true;
        }

        if (!opts.unclassified_output_filename.is_empty()) {
            if (opts.paired_end_processing) {
                let fields = SplitString(&opts.unclassified_output_filename, "#");
                if (fields.len() != 2) {
                    return Err("Paired filename format missing or extra '#'".into());
                }
                outputs.unclassified_output1 = Some(BufWriter::new(File::create(format!(
                    "{}_1{}",
                    fields[0], fields[1]
                ))?));
                outputs.unclassified_output2 = Some(BufWriter::new(File::create(format!(
                    "{}_2{}",
                    fields[0], fields[1]
                ))?));
            } else {
                outputs.unclassified_output1 = Some(BufWriter::new(File::create(
                    &opts.unclassified_output_filename,
                )?));
            }
            outputs.printing_sequences = true;
        }

        if (!opts.kraken_output_filename.is_empty() && opts.kraken_output_filename != "-") {
            outputs.kraken_output =
                Some(BufWriter::new(File::create(&opts.kraken_output_filename)?));
        }

        outputs.initialized = true;
    }

    Ok(())
}

fn add_hitlist_string(oss: &mut String, taxa: &[taxid_t], taxonomy: &Taxonomy) {
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in &taxa[1..] {
        if (code == last_code) {
            code_count += 1;
        } else {
            if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON) {
                if (last_code == AMBIGUOUS_SPAN_TAXON) {
                    oss.push_str(&format!("A:{} ", code_count));
                } else {
                    let ext_code = taxonomy.nodes()[last_code as usize].external_id;
                    oss.push_str(&format!("{}:{} ", ext_code, code_count));
                }
            } else {
                // mate pair / reading frame marker
                oss.push_str(if (last_code == MATE_PAIR_BORDER_TAXON) {
                    "|:| "
                } else {
                    "-:- "
                });
            }
            code_count = 1;
            last_code = code;
        }
    }
    // final block
    if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON) {
        if (last_code == AMBIGUOUS_SPAN_TAXON) {
            oss.push_str(&format!("A:{}", code_count));
        } else {
            let ext_code = taxonomy.nodes()[last_code as usize].external_id;
            oss.push_str(&format!("{}:{}", ext_code, code_count));
        }
    } else {
        oss.push_str(if (last_code == MATE_PAIR_BORDER_TAXON) {
            "|:|"
        } else {
            "-:-"
        });
    }
}

fn resolve_tree(
    hit_counts: &TaxonCountsT,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> taxid_t {
    let required_score = (opts.confidence_threshold * (total_minimizers as f64)).ceil() as u32;

    let mut max_taxon = 0;
    let mut max_score = 0;

    // sum each taxon's LTR path
    for &taxon in hit_counts.keys() {
        let mut score = 0;
        for (&taxon2, &cnt2) in hit_counts.iter() {
            if (taxonomy.is_ancestor_of(taxon2, taxon)) {
                score += cnt2;
            }
        }
        if (score > max_score) {
            max_score = score;
            max_taxon = taxon;
        } else if (score == max_score) {
            max_taxon = taxonomy.lowest_common_ancestor(max_taxon, taxon);
        }
    }

    // reset max_score
    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);
    let mut current_taxon = max_taxon;

    while (current_taxon != 0 && max_score < required_score) {
        // sum hits in that clade
        let mut clade_score = 0;
        for (&t, &c) in hit_counts.iter() {
            if (taxonomy.is_ancestor_of(current_taxon, t)) {
                clade_score += c;
            }
        }
        if (clade_score >= required_score) {
            return current_taxon;
        } else {
            current_taxon = taxonomy.nodes()[current_taxon as usize].parent_id;
        }
    }

    current_taxon
}

pub fn classify_sequence(
    dna: &mut Sequence,
    dna2: &mut Sequence,
    koss: &mut String,
    hash: &KeyValueStore,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<taxid_t>,
    hit_counts: &mut HashMap<taxid_t, u32>,
    tx_frames: &mut Vec<String>,
    curr_taxon_counts: &mut TaxonCountsT,
) -> taxid_t {
    taxa.clear();
    hit_counts.clear();

    let mut call: taxid_t = 0;
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups: i64 = 0;

    // Possibly mask low-quality bases
    if opts.minimum_quality_score > 0 {
        mask_low_quality_bases(dna, opts.minimum_quality_score);
        if opts.paired_end_processing {
            mask_low_quality_bases(dna2, opts.minimum_quality_score);
        }
    }

    // Translated search pre-processing
    if opts.use_translated_search {
        TranslateToAllFrames(&dna.seq, tx_frames);
    }

    'outer: for mate_num in 0..2 {
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        for frame_idx in 0..frame_ct {
            if opts.use_translated_search {
                scanner.load_sequence(&tx_frames[frame_idx]);
            } else {
                let seq = if mate_num == 0 { &dna.seq } else { &dna2.seq };
                scanner.load_sequence(seq);
            }

            let mut last_minimizer = u64::MAX;
            let mut last_taxon = u64::MAX;
            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon = if scanner.is_ambiguous() {
                    AMBIGUOUS_SPAN_TAXON
                } else {
                    if *minimizer != last_minimizer {
                        let mut found_taxon: taxid_t = 0;
                        if idx_opts.minimum_acceptable_hash_value == 0
                            || murmur_hash3(*minimizer) >= idx_opts.minimum_acceptable_hash_value
                        {
                            found_taxon = hash.get(*minimizer);
                        }
                        last_taxon = found_taxon;
                        last_minimizer = *minimizer;
                        if found_taxon != 0 {
                            minimizer_hit_groups += 1;
                            // Register kmer in curr_taxon_counts
                            curr_taxon_counts
                                .entry(found_taxon)
                                .or_insert_with(TaxonCounter::default)
                                .add_kmer(scanner.last_minimizer());
                        }
                        found_taxon
                    } else {
                        last_taxon
                    }
                };

                if taxon != 0 {
                    if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                        call = taxon;
                        break 'outer;
                    }
                    *hit_counts.entry(taxon).or_insert(0) += 1;
                }
                taxa.push(taxon);
            }

            if opts.use_translated_search && frame_idx != frame_ct - 1 {
                taxa.push(READING_FRAME_BORDER_TAXON);
            }
        }

        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    // Compute call if not quick-mode decided
    if call == 0 {
        let mut total_kmers = taxa.len();
        if opts.paired_end_processing {
            total_kmers -= 1; // account for mate pair marker
        }
        if opts.use_translated_search {
            // account for reading frame markers (2 if single-end, 4 if paired)
            total_kmers -= if opts.paired_end_processing { 4 } else { 2 };
        }

        call = ResolveTree(hit_counts, taxonomy, total_kmers, opts);
        if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
            call = 0;
        }
    }

    // Update stats and classification counters
    if call != 0 {
        stats.total_classified += 1;
        curr_taxon_counts
            .entry(call)
            .or_insert_with(TaxonCounter::default)
            .increment_read_count();
    }

    // Print result line
    if call != 0 {
        write!(koss, "C\t").unwrap();
    } else {
        write!(koss, "U\t").unwrap();
    }

    let read_id = if opts.paired_end_processing {
        TrimPairInfo(&dna.id)
    } else {
        dna.id.clone()
    };
    write!(koss, "{}\t", read_id).unwrap();

    let ext_call = taxonomy.nodes()[call].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            taxonomy.name_data_at_offset(taxonomy.nodes()[call].name_offset)
        } else {
            "unclassified"
        };
        write!(koss, "{} (taxid {})\t", name, ext_call).unwrap();
    } else {
        write!(koss, "{}\t", ext_call).unwrap();
    }

    if !opts.paired_end_processing {
        write!(koss, "{}\t", dna.seq.len()).unwrap();
    } else {
        write!(koss, "{}|{}\t", dna.seq.len(), dna2.seq.len()).unwrap();
    }

    if opts.quick_mode {
        write!(koss, "{}:Q\n", ext_call).unwrap();
    } else {
        if taxa.is_empty() {
            write!(koss, "0:0\n").unwrap();
        } else {
            AddHitlistString(koss, taxa, taxonomy);
            writeln!(koss).unwrap();
        }
    }

    call
}

fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if (sz <= 2) {
        return id.to_string();
    }
    let bytes = id.as_bytes();
    if (bytes[sz - 2] == b'/' && (bytes[sz - 1] == b'1' || bytes[sz - 1] == b'2')) {
        id[..sz - 2].to_string()
    } else {
        id.to_string()
    }
}

fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: &mut KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    total_taxon_counters: &mut taxon_counters_t,
) {
    // Handle file openings
    let fptr1: Box<dyn Read> = match filename1 {
        Some(name) => Box::new(File::open(name).expect("Failed to open file")),
        None => Box::new(std::io::stdin()),
    };

    let fptr2: Option<Box<dyn Read>> = if opts.paired_end_processing && !opts.single_file_pairs {
        filename2
            .map(|name| Box::new(File::open(name).expect("Failed to open file")) as Box<dyn Read>)
    } else {
        None
    };

    // Implement the priority queue
    struct OutputData {
        block_id: u64,
        kraken_str: String,
        classified_out1_str: String,
        classified_out2_str: String,
        unclassified_out1_str: String,
        unclassified_out2_str: String,
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

    impl Eq for OutputData {}

    impl PartialEq for OutputData {
        fn eq(&self, other: &Self) -> bool {
            self.block_id == other.block_id
        }
    }

    let output_queue: Arc<Mutex<BinaryHeap<OutputData>>> = Arc::new(Mutex::new(BinaryHeap::new()));
    let next_input_block_id: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));
    let next_output_block_id: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));

    // Handle concurrency with threads
    let num_threads = opts.num_threads;
    let mut threads: Vec<std::thread::JoinHandle<()>> = Vec::new();

    for _ in 0..num_threads {
        let hash_clone = hash.clone(); // Assuming KeyValueStore implements Clone
        let tax_clone = tax.clone(); // Assuming Taxonomy implements Clone
        let idx_opts_clone = idx_opts.clone(); // Assuming IndexOptions implements Clone
        let opts_clone = opts.clone(); // Assuming Options implements Clone
        let stats_clone = stats.clone(); // Assuming ClassificationStats implements Clone
        let outputs_clone = outputs.clone(); // Assuming OutputStreamData implements Clone
        let total_taxon_counters_clone = total_taxon_counters.clone(); // Assuming taxon_counters_t implements Clone
        let output_queue_clone = Arc::clone(&output_queue);
        let next_input_block_id_clone = Arc::clone(&next_input_block_id);

        threads.push(thread::spawn(move || {
            // Thread-local data
            let mut scanner = MinimizerScanner::new(
                idx_opts_clone.k,
                idx_opts_clone.l,
                idx_opts_clone.spaced_seed_mask,
                idx_opts_clone.dna_db,
                idx_opts_clone.toggle_mask,
                idx_opts_clone.revcom_version,
            );
            let mut taxa = Vec::new();
            let mut hit_counts = taxon_counts_t::new();
            let mut kraken_oss = String::new();
            let mut c1_oss = String::new();
            let mut c2_oss = String::new();
            let mut u1_oss = String::new();
            let mut u2_oss = String::new();
            let mut translated_frames = vec![String::new(); 6];
            let mut thread_taxon_counters = taxon_counters_t::new();

            loop {
                let block_id = {
                    let mut lock = next_input_block_id_clone.lock().unwrap();
                    let id = *lock;
                    *lock += 1;
                    id
                };

                // Load block of sequences
                let mut reader1 = BatchSequenceReader::new();
                let mut reader2 = BatchSequenceReader::new();
                let ok_read = if !opts.paired_end_processing {
                    reader1.load_block(&*fptr1, 3 * 1024 * 1024)
                } else if !opts.single_file_pairs {
                    reader1.load_batch(&*fptr1, NUM_FRAGMENTS_PER_THREAD);
                    reader2.load_batch(fptr2.as_ref().unwrap(), NUM_FRAGMENTS_PER_THREAD)
                } else {
                    reader1.load_batch(&*fptr1, NUM_FRAGMENTS_PER_THREAD * 2);
                };

                if !ok_read {
                    break;
                }

                // Process sequences
                while reader1.next_sequence(&mut seq1) {
                    if opts.paired_end_processing {
                        reader2.next_sequence(&mut seq2);
                    }

                    // Process sequence classification
                    let call = classify_sequence(
                        &seq1,
                        if opts.paired_end_processing {
                            &seq2
                        } else {
                            &mut Sequence::default()
                        },
                        &mut kraken_oss,
                        &hash_clone,
                        &tax_clone,
                        &idx_opts_clone,
                        &opts_clone,
                        &mut stats_clone,
                        &mut scanner,
                        &mut taxa,
                        &mut hit_counts,
                        &mut translated_frames,
                        &mut thread_taxon_counters,
                    );

                    // Update output strings based on call
                    if call != 0 {
                        // Classified
                        seq1.header.push_str(&format!(
                            " kraken:taxid|{}",
                            tax.nodes()[call as usize].external_id
                        ));
                        c1_oss.push_str(&seq1.to_string());
                        if opts.paired_end_processing {
                            seq2.header.push_str(&format!(
                                " kraken:taxid|{}",
                                tax.nodes()[call as usize].external_id
                            ));
                            c2_oss.push_str(&seq2.to_string());
                        }
                    } else {
                        // Unclassified
                        u1_oss.push_str(&seq1.to_string());
                        if opts.paired_end_processing {
                            u2_oss.push_str(&seq2.to_string());
                        }
                    }

                    // Update statistics
                    stats_clone.total_sequences += 1;
                    stats_clone.total_bases += seq1.seq.len();
                    if opts.paired_end_processing {
                        stats_clone.total_bases += seq2.seq.len();
                    }
                    if call != 0 {
                        stats_clone.total_classified += 1;
                    }
                }

                // Push output data to queue
                let mut out_data = OutputData {
                    block_id,
                    kraken_str: kraken_oss.clone(),
                    classified_out1_str: c1_oss.clone(),
                    classified_out2_str: c2_oss.clone(),
                    unclassified_out1_str: u1_oss.clone(),
                    unclassified_out2_str: u2_oss.clone(),
                };
                output_queue_clone.lock().unwrap().push(out_data);

                // Reset thread-local data for next block
                kraken_oss.clear();
                c1_oss.clear();
                c2_oss.clear();
                u1_oss.clear();
                u2_oss.clear();
            }
        }));
    }

    // Join all threads
    for handle in threads {
        handle.join().unwrap();
    }

    // Output data in order
    let mut output_queue = output_queue.lock().unwrap();
    while let Some(out_data) = output_queue.pop() {
        if out_data.block_id == *next_output_block_id.lock().unwrap() {
            // Write to output streams
            if let Some(ref mut kraken_output) = outputs.kraken_output {
                kraken_output
                    .write_all(out_data.kraken_str.as_bytes())
                    .unwrap();
            }
            if let Some(ref mut classified_output1) = outputs.classified_output1 {
                classified_output1
                    .write_all(out_data.classified_out1_str.as_bytes())
                    .unwrap();
            }
            // Similarly write other output strings
            *next_output_block_id.lock().unwrap() += 1;
        } else {
            // Need to handle out of order data
            // For now, push back if not in order
            output_queue.push(out_data);
            break;
        }
    }
}

fn report_stats(start: Instant, stats: &ClassificationStats) {
    let elapsed = start.elapsed();
    let seconds = elapsed.as_secs_f64();
    let total_unclassified = stats.total_sequences - stats.total_classified;

    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1.0e6,
        seconds,
        stats.total_sequences as f64 / 1.0e3 / (seconds / 60.0),
        stats.total_bases as f64 / 1.0e6 / (seconds / 60.0)
    );
    eprintln!(
        "  {} sequences classified ({:.2}%)",
        stats.total_classified,
        stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
    );
    eprintln!(
        "  {} sequences unclassified ({:.2}%)",
        total_unclassified,
        total_unclassified as f64 * 100.0 / stats.total_sequences as f64
    );
}

fn parse_command_line() -> Result<Options, Box<dyn Error>> {
    let matches = App::new("classify")
        .arg(
            Arg::with_name("H")
                .short("H")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("t")
                .short("t")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("o")
                .short("o")
                .takes_value(true)
                .required(true),
        )
        // ... add other args similarly ...
        .get_matches();

    let mut opts = Options::default();
    opts.index_filename = matches.value_of("H").unwrap().to_string();
    opts.taxonomy_filename = matches.value_of("t").unwrap().to_string();
    opts.options_filename = matches.value_of("o").unwrap().to_string();
    // Parse other options similarly
    Ok(opts)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut opts = parse_command_line()?;

    // Load index options
    // load idx_opts from file
    let idx_opts = IndexOptions {
        k: 31,
        l: 15,
        spaced_seed_mask: 0,
        dna_db: true,
        toggle_mask: 0,
        revcom_version: 2,
        minimum_acceptable_hash_value: 0,
    };

    // load taxonomy
    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping)?;

    // load hash
    let hash_ptr =
        compact_hash::CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping)?;

    let mut stats = ClassificationStats::default();
    let mut outputs = OutputStreamData::default();
    let mut taxon_counters = TaxonCountersT::new();

    let start = Instant::now();

    // Process input files
    // If no input files, process stdin
    let args: Vec<String> = std::env::args().collect();
    // This is simplified logic from original code
    let input_files = &args[1..];

    if (input_files.is_empty()) {
        if (opts.paired_end_processing && !opts.single_file_pairs) {
            return Err("paired end processing used with no files specified".into());
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
        while (i < input_files.len()) {
            if (opts.paired_end_processing && !opts.single_file_pairs) {
                if (i + 1 == input_files.len()) {
                    return Err("paired end processing used with unpaired file".into());
                }
                process_files(
                    Some(&input_files[i]),
                    Some(&input_files[i + 1]),
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
                    Some(&input_files[i]),
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

    report_stats(start, &stats);

    if (!opts.report_filename.is_empty()) {
        let total_unclassified = stats.total_sequences - stats.total_classified;
        if (opts.mpa_style_report) {
            report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &taxon_counters,
            )?;
        } else {
            ReportKrakenStyle(
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

fn classify_reads(reads: Vec<String>, db: Arc<Mutex<HashMap<String, String>>>) -> Vec<String> {
    let mut results = Vec::new();
    let mut handles = vec![];

    for read in reads {
        let db = Arc::clone(&db);
        let handle = thread::spawn(move || {
            let db = db.lock().unwrap();
            // ...existing code...
            results.push(result);
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    results
}
