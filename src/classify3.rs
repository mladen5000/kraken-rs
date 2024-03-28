/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::{Arc, Barrier, Mutex};
use std::time::{Duration, Instant, SystemTime};

use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{IndexOptions, TaxonCounts};
use crate::mmscanner::MinimizerScanner;
use crate::reports::TaxonCounters;
use crate::seqreader::{BatchSequenceReader, Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy; // Add this import statement

type TaxId = u64;
const TAXID_MAX: TaxId = std::u64::MAX;

const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
const MATE_PAIR_BORDER_TAXON: TaxId = TaxId::MAX;
const READING_FRAME_BORDER_TAXON: TaxId = TaxId::MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxId = TaxId::MAX - 2;

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
    minimum_quality_score: u8,
    minimum_hit_groups: usize,
    use_memory_mapping: bool,
    match_input_order: bool,
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

struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}

fn load_index_options(options_filename: &str) -> IndexOptions {
    let mut idx_opts = IndexOptions {
        k: 0,
        l: 0,
        spaced_seed_mask: 0,
        dna_db: false,
        minimum_acceptable_hash_value: 0, // Fix the type mismatch by assigning `None` as a `u64` value
        toggle_mask: todo!(),
        revcom_version: todo!(),
        db_version: todo!(),
        db_type: todo!(),
    };

    let file = File::open(options_filename).expect("Unable to open options file");
    let mut reader = BufReader::new(file);
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .expect("Unable to read options file");

    let opts_filesize = buffer.len();
    idx_opts.k = match std::mem::size_of::<usize>() {
        4 => u32::from_le_bytes([buffer[0], buffer[1], buffer[2], buffer[3]]) as usize,
        8 => u64::from_le_bytes([
            buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6], buffer[7],
        ]) as usize,
        _ => panic!("Unsupported architecture"),
    };

    idx_opts.l = match std::mem::size_of::<usize>() {
        4 => u32::from_le_bytes([buffer[8], buffer[9], buffer[10], buffer[11]]) as usize,
        8 => u64::from_le_bytes([
            buffer[12], buffer[13], buffer[14], buffer[15], buffer[16], buffer[17], buffer[18],
            buffer[19],
        ]) as usize,
        _ => panic!("Unsupported architecture"),
    };
    idx_opts.spaced_seed_mask = u64::from_le_bytes([
        buffer[8], buffer[9], buffer[10], buffer[11], buffer[12], buffer[13], buffer[14],
        buffer[15],
    ]);
    idx_opts.dna_db = buffer[16] != 0;
    if opts_filesize >= 25 {
        idx_opts.minimum_acceptable_hash_value = u64::from_le_bytes([
            buffer[17], buffer[18], buffer[19], buffer[20], buffer[21], buffer[22], buffer[23],
            buffer[24],
        ]);
    }

    idx_opts
}

fn process_input_files(
    opts: &Options,
    idx_opts: &IndexOptions,
    taxonomy: &Taxonomy,
    hash_ptr: &CompactHashTable,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
) {
    let args: Vec<String> = std::env::args().collect();

    if args.len() == 1 {
        if opts.paired_end_processing && !opts.single_file_pairs {
            panic!("Paired end processing used with no files specified");
        }
        process_files(
            None,
            None,
            hash_ptr,
            taxonomy,
            idx_opts,
            opts,
            stats,
            outputs,
            taxon_counters,
        );
    } else {
        for mut i in 1..args.len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == args.len() {
                    panic!("Paired end processing used with unpaired file");
                }
                process_files(
                    Some(&args[i]),
                    Some(&args[i + 1]),
                    hash_ptr,
                    taxonomy,
                    idx_opts,
                    opts,
                    stats,
                    outputs,
                    taxon_counters,
                );
                i += 1;
            } else {
                process_files(
                    Some(&args[i]),
                    None,
                    hash_ptr,
                    taxonomy,
                    idx_opts,
                    opts,
                    stats,
                    outputs,
                    taxon_counters,
                );
            }
        }
    }
}

fn report_stats(time1: Duration, time2: Duration, stats: &ClassificationStats) {
    let time_diff = time2 - time1;
    let seconds = time_diff.as_secs_f64();

    let total_unclassified = stats.total_sequences - stats.total_classified;

    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }

    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1.0e6,
        seconds,
        stats.total_sequences as f64 / 1.0e3 / (seconds / 60.0),
        stats.total_bases as f64 / 1.0e6 / (seconds / 60.0)
    );

    eprintln!(
        " {} sequences classified ({:.2}%)",
        stats.total_classified,
        stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
    );

    eprintln!(
        " {} sequences unclassified ({:.2}%)",
        total_unclassified,
        total_unclassified as f64 * 100.0 / stats.total_sequences as f64
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
    total_taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
) {
    let (fptr1, fptr2) = open_input_files(filename1, filename2, opts);

    let mut output_queue = initialize_output_queue();
    let mut next_input_block_id = 0;
    let mut next_output_block_id = 0;
    let output_lock = Arc::new(Mutex::new(()));

    let thread_count = get_thread_count();
    let barrier = Arc::new(Barrier::new(thread_count));

    (0..thread_count).into_par_iter().for_each(|_| {
        let mut scanner = MinimizerScanner::new(
            idx.opts.k,
            idx_opts.l,
            &idx_opts.spaced_seed_mask,
            &idx_opts.dna_db,
            &idx_opts.toggle_mask,
            idx_opts.revcom_version,

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
        let mut reader1 = BatchSequenceReader::new();
        let mut reader2 = BatchSequenceReader::new();
        let mut seq1 = Sequence::new();
        let mut seq2 = Sequence::new();
        let mut thread_taxon_counters = HashMap::new();

        loop {
            thread_stats.reset();

            let (ok_read, block_id) = read_input_block(
                &mut reader1,
                &mut reader2,
                &mut next_input_block_id,
                &fptr1,
                &fptr2,
                opts,
            );

            if !ok_read {
                break;
            }

            reset_output_strings(
                &mut kraken_oss,
                &mut c1_oss,
                &mut c2_oss,
                &mut u1_oss,
                &mut u2_oss,
                &mut thread_taxon_counters,
            );

            while let Some((seq1, seq2)) =
                read_sequence_pair(&mut reader1, &mut reader2, &mut seq1, &mut seq2, opts)
            {
                process_sequence_pair(
                    &mut seq1,
                    &mut seq2,
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
                    &mut kraken_oss,
                    &mut c1_oss,
                    &mut c2_oss,
                    &mut u1_oss,
                    &mut u2_oss,
                );
            }

            update_global_stats(stats, &thread_stats);
            print_progress(stats);

            initialize_outputs_if_needed(opts, outputs, reader1.file_format());

            let out_data =
                prepare_output_data(block_id, &kraken_oss, &c1_oss, &c2_oss, &u1_oss, &u2_oss);

            output_queue.lock().unwrap().push(out_data);

            update_total_taxon_counters(total_taxon_counters, &mut thread_taxon_counters);

            process_output_queue(
                &output_queue,
                &mut next_output_block_id,
                outputs,
                &output_lock,
            );
        }
    });

    close_input_files(fptr1, fptr2, opts);
    flush_output_streams(outputs);
}

use std::cmp::Reverse;

fn open_input_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    opts: &Options,
) -> (Option<Box<dyn Read>>, Option<Box<dyn Read>>) {
    let fptr1 = match filename1 {
        Some(file) => Some(Box::new(BufReader::new(File::open(file).unwrap())) as Box<dyn Read>),
        None => None,
    };

    let fptr2 = if opts.paired_end_processing && !opts.single_file_pairs {
        match filename2 {
            Some(file) => {
                Some(Box::new(BufReader::new(File::open(file).unwrap())) as Box<dyn Read>)
            }
            None => panic!("Second input file is required for paired-end processing."),
        }
    } else {
        None
    };

    (fptr1, fptr2)
}

fn initialize_output_queue() -> Arc<Mutex<BinaryHeap<Reverse<OutputData>>>> {
    let comparator = |a: &OutputData, b: &OutputData| a.block_id.cmp(&b.block_id);
    Arc::new(Mutex::new(BinaryHeap::with_capacity_by(1024, comparator)))
}

fn get_thread_count() -> usize {
    num_cpus::get()
}

fn read_input_block(
    reader1: &mut BatchSequenceReader,
    reader2: &mut BatchSequenceReader,
    next_input_block_id: &mut u64,
    fptr1: &Option<Box<dyn Read>>,
    fptr2: &Option<Box<dyn Read>>,
    opts: &Options,
) -> (bool, u64) {
    let mut ok_read = false;
    let block_id = *next_input_block_id;

    if let Some(ref mut fptr1) = fptr1 {
        if !opts.paired_end_processing {
            // Unpaired data? Just read in a sized block
            ok_read = reader1.load_block(fptr1, 3 * 1024 * 1024);
        } else if !opts.single_file_pairs {
            // Paired data in 2 files? Read a line-counted batch from each file.
            ok_read = reader1.load_batch(fptr1, NUM_FRAGMENTS_PER_THREAD);
            if ok_read && opts.paired_end_processing {
                if let Some(ref mut fptr2) = fptr2 {
                    ok_read = reader2.load_batch(fptr2, NUM_FRAGMENTS_PER_THREAD);
                } else {
                    panic!("Second input file is required for paired-end processing.");
                }
            }
        } else {
            let mut num_fragments = NUM_FRAGMENTS_PER_THREAD * 2;
            // Ensure fragment count is even - just in case the above line is changed
            if num_fragments % 2 == 1 {
                num_fragments += 1;
            }
            ok_read = reader1.load_batch(fptr1, num_fragments);
        }
    }

    *next_input_block_id += 1;
    (ok_read, block_id)
}

fn reset_output_strings(
    kraken_oss: &mut String,
    c1_oss: &mut String,
    c2_oss: &mut String,
    u1_oss: &mut String,
    u2_oss: &mut String,
    thread_taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
) {
    kraken_oss.clear();
    c1_oss.clear();
    c2_oss.clear();
    u1_oss.clear();
    u2_oss.clear();
    thread_taxon_counters.clear();
}

fn read_sequence_pair(
    reader1: &mut BatchSequenceReader,
    reader2: &mut BatchSequenceReader,
    seq1: &mut Sequence,
    seq2: &mut Sequence,
    opts: &Options,
) -> Option<(Sequence, Sequence)> {
    let valid_fragment = reader1.next_sequence();
    if !valid_fragment {
        return None;
    }

    if opts.paired_end_processing {
        let valid_fragment = if opts.single_file_pairs {
            reader1.next_sequence()
        } else {
            reader2.next_sequence()
        };

        if !valid_fragment {
            return None;
        }
    }

    Some((seq1.clone(), seq2.clone()))
}

fn process_sequence_pair(
    seq1: &mut Sequence,
    seq2: &mut Sequence,
    hash: &dyn KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    thread_stats: &mut ClassificationStats,
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<TaxId>,
    hit_counts: &mut HashMap<TaxId, u32>,
    translated_frames: &mut [String],
    thread_taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
    kraken_oss: &mut String,
    c1_oss: &mut String,
    c2_oss: &mut String,
    u1_oss: &mut String,
    u2_oss: &mut String,
) {
    thread_stats.total_sequences += 1;

    if opts.minimum_quality_score > 0 {
        mask_low_quality_bases(seq1, opts.minimum_quality_score);
        if opts.paired_end_processing {
            mask_low_quality_bases(seq2, opts.minimum_quality_score);
        }
    }

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
    );

    if call != 0 {
        let buffer = format!(" kraken:taxid|{}", tax.nodes()[call as usize].external_id);
        seq1.header += &buffer;
        c1_oss.push_str(&seq1.to_string());
        if opts.paired_end_processing {
            seq2.header += &buffer;
            c2_oss.push_str(&seq2.to_string());
        }
    } else {
        u1_oss.push_str(&seq1.to_string());
        if opts.paired_end_processing {
            u2_oss.push_str(&seq2.to_string());
        }
    }

    thread_stats.total_bases += seq1.seq.len();
    if opts.paired_end_processing {
        thread_stats.total_bases += seq2.seq.len();
    }
}

fn update_global_stats(stats: &mut ClassificationStats, thread_stats: &ClassificationStats) {
    stats.total_sequences += thread_stats.total_sequences;
    stats.total_bases += thread_stats.total_bases;
    stats.total_classified += thread_stats.total_classified;
}

fn print_progress(stats: &ClassificationStats) {
    if atty::is(atty::Stream::Stderr) {
        eprint!(
            "\rProcessed {} sequences ({} bp) ...",
            stats.total_sequences, stats.total_bases
        );
    }
}

fn initialize_outputs_if_needed(
    opts: &Options,
    outputs: &mut OutputStreamData,
    file_format: SequenceFormat,
) {
    if !outputs.initialized {
        initialize_outputs(opts, outputs, file_format);
    }
}

fn prepare_output_data(
    block_id: u64,
    kraken_oss: &str,
    c1_oss: &str,
    c2_oss: &str,
    u1_oss: &str,
    u2_oss: &str,
) -> OutputData {
    OutputData {
        block_id,
        kraken_str: kraken_oss.to_string(),
        classified_out1_str: c1_oss.to_string(),
        classified_out2_str: c2_oss.to_string(),
        unclassified_out1_str: u1_oss.to_string(),
        unclassified_out2_str: u2_oss.to_string(),
    }
}

fn update_total_taxon_counters(
    total_taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
    thread_taxon_counters: &mut HashMap<TaxId, TaxonCounters>,
) {
    for (taxid, counter) in thread_taxon_counters.drain() {
        total_taxon_counters
            .entry(taxid)
            .or_default()
            .merge(counter);
    }
}

fn process_output_queue(
    output_queue: &Arc<Mutex<BinaryHeap<Reverse<OutputData>>>>,
    next_output_block_id: &mut u64,
    outputs: &mut OutputStreamData,
    output_lock: &Arc<Mutex<()>>,
) {
    loop {
        let mut output_data = None;
        {
            let mut queue = output_queue.lock().unwrap();
            if let Some(data) = queue.peek() {
                if data.0.block_id == *next_output_block_id {
                    output_data = Some(queue.pop().unwrap().0);
                    *next_output_block_id += 1;
                }
            }
        }

        if let Some(data) = output_data {
            let _lock = output_lock.lock().unwrap();
            if let Some(ref mut kraken_output) = outputs.kraken_output {
                kraken_output.write_all(data.kraken_str.as_bytes()).unwrap();
            }
            if let Some(ref mut classified_output1) = outputs.classified_output1 {
                classified_output1
                    .write_all(data.classified_out1_str.as_bytes())
                    .unwrap();
            }
            if let Some(ref mut classified_output2) = outputs.classified_output2 {
                classified_output2
                    .write_all(data.classified_out2_str.as_bytes())
                    .unwrap();
            }
            if let Some(ref mut unclassified_output1) = outputs.unclassified_output1 {
                unclassified_output1
                    .write_all(data.unclassified_out1_str.as_bytes())
                    .unwrap();
            }
            if let Some(ref mut unclassified_output2) = outputs.unclassified_output2 {
                unclassified_output2
                    .write_all(data.unclassified_out2_str.as_bytes())
                    .unwrap();
            }
        } else {
            break;
        }
    }
}

fn close_input_files(fptr1: Option<Box<dyn Read>>, fptr2: Option<Box<dyn Read>>, _opts: &Options) {
    if let Some(fptr) = fptr1 {
        drop(fptr);
    }
    if let Some(fptr) = fptr2 {
        drop(fptr);
    }
}

fn flush_output_streams(outputs: &mut OutputStreamData) {
    if let Some(ref mut kraken_output) = outputs.kraken_output {
        kraken_output.flush().unwrap();
    }
    if let Some(ref mut classified_output1) = outputs.classified_output1 {
        classified_output1.flush().unwrap();
    }
    if let Some(ref mut classified_output2) = outputs.classified_output2 {
        classified_output2.flush().unwrap();
    }
    if let Some(ref mut unclassified_output1) = outputs.unclassified_output1 {
        unclassified_output1.flush().unwrap();
    }
    if let Some(ref mut unclassified_output2) = outputs.unclassified_output2 {
        unclassified_output2.flush().unwrap();
    }
}

fn resolve_tree(
    hit_counts: &mut HashMap<TaxId, u32>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (opts.confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for (&taxon, _) in hit_counts.iter() {
        let mut score = 0;
        for (&taxon2, &count) in hit_counts.iter() {
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
    max_score = hit_counts[&max_taxon];

    // We probably have a call w/o required support (unless LCA resolved tie)
    while max_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts.iter() {
            // Add to score if taxon in max_taxon's clade
            if taxonomy.is_a_ancestor_of_b(max_taxon, taxon) {
                max_score += count;
            }
        }

        // Score is now sum of hits at max_taxon and w/in max_taxon clade
        if max_score >= required_score {
            // Kill loop and return, we've got enough support here
            return max_taxon;
        } else {
            // Run up tree until confidence threshold is met
            // Run off tree if required score isn't met
            max_taxon = taxonomy.nodes[&max_taxon].parent_id;
        }
    }

    max_taxon
}

/// Consistent translation
fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        id.to_string()
    } else if &id[sz - 2..] == "/1" || &id[sz - 2..] == "/2" {
        id[..sz - 2].to_string()
    } else {
        id.to_string()
    }
}

use crate::kv_store::{murmurhash3, KeyValueStore};

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
    hit_counts: &mut HashMap<TaxId, u32>,
    tx_frames: &mut Vec<String>,
    curr_taxon_counts: &mut HashMap<TaxId, TaxonCounters>,
) -> TaxId {
    let mut call = 0;
    taxa.clear();
    hit_counts.clear();
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups = 0;

    for mate_num in 0..2 {
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        if opts.use_translated_search {
            translate_to_all_frames(if mate_num == 0 { &dna.seq } else { &dna2.seq });
        }

        for frame_idx in 0..frame_ct {
            scan_minimizers(
                opts,
                idx_opts,
                scanner,
                hash,
                taxa,
                hit_counts,
                &mut minimizer_hit_groups,
                tx_frames,
                mate_num,
                frame_idx,
                dna,
                dna2,
                curr_taxon_counts,
                &mut call,
            );

            if opts.use_translated_search && frame_idx != 5 {
                taxa.push(READING_FRAME_BORDER_TAXON);
            }
        }

        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    let total_kmers = calculate_total_kmers(opts, taxa);
    call = resolve_tree(hit_counts, taxonomy, total_kmers, opts);
    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    update_stats(call, stats, curr_taxon_counts);
    write_classification_output(call, opts, taxonomy, dna, dna2, koss, taxa);

    call
}

fn scan_minimizers(
    opts: &Options,
    idx_opts: &IndexOptions,
    scanner: &mut MinimizerScanner,
    hash: &dyn KeyValueStore,
    taxa: &mut Vec<TaxId>,
    hit_counts: &mut HashMap<TaxId, u32>,
    minimizer_hit_groups: &mut i64,
    tx_frames: &[String],
    mate_num: usize,
    frame_idx: usize,
    dna: &Sequence,
    dna2: &Sequence,
    curr_taxon_counts: &mut HashMap<TaxId, TaxonCounters>,
    call: &mut TaxId,
) {
    if opts.use_translated_search {
        // DEbUG: I added start + finish parameters of 0 and 100
        scanner.load_sequence(&tx_frames[frame_idx], 0, 100);
    } else {
        scanner.load_sequence(if mate_num == 0 { &dna.seq } else { &dna2.seq }, 0, 100);
    }

    let mut last_minimizer = u64::MAX;
    let mut last_taxon = TAXID_MAX;

    while let Some(minimizer_ptr) = scanner.next_minimizer() {
        let taxon = process_minimizer(
            scanner,
            idx_opts,
            minimizer_ptr,
            hash,
            &mut last_minimizer,
            &mut last_taxon,
            minimizer_hit_groups,
            curr_taxon_counts,
        );

        if taxon != 0 {
            if opts.quick_mode && *minimizer_hit_groups >= opts.minimum_hit_groups {
                *call = taxon;
                return;
            }
            *hit_counts.entry(taxon).or_default() += 1;
        }
        taxa.push(taxon);
    }
}

fn process_minimizer(
    scanner: &MinimizerScanner,
    idx_opts: &IndexOptions,
    minimizer_ptr: &u64,
    hash: &dyn KeyValueStore,
    last_minimizer: &mut u64,
    last_taxon: &mut TaxId,
    minimizer_hit_groups: &mut i64,
    curr_taxon_counts: &mut HashMap<TaxId, TaxonCounters>,
) -> TaxId {
    if scanner.is_ambiguous() {
        return AMBIGUOUS_SPAN_TAXON;
    }

    let mut taxon = 0;
    if *minimizer_ptr != *last_minimizer {
        let mut skip_lookup = false;
        if let Some(min_hash_value) = idx_opts.minimum_acceptable_hash_value {
            if murmur_hash3(*minimizer_ptr) < min_hash_value {
                skip_lookup = true;
            }
        }
        if !skip_lookup {
            taxon = hash.get(*minimizer_ptr);
        }
        *last_taxon = taxon;
        *last_minimizer = *minimizer_ptr;
        if taxon != 0 {
            *minimizer_hit_groups += 1;
            curr_taxon_counts
                .entry(taxon)
                .or_default()
                .add_kmer(scanner.last_minimizer());
        }
    } else {
        taxon = *last_taxon;
    }

    taxon
}

fn calculate_total_kmers(opts: &Options, taxa: &[TaxId]) -> usize {
    let mut total_kmers = taxa.len();
    if opts.paired_end_processing {
        total_kmers -= 1;
    }
    if opts.use_translated_search {
        total_kmers -= if opts.paired_end_processing { 4 } else { 2 };
    }
    total_kmers
}

fn update_stats(
    call: TaxId,
    stats: &mut ClassificationStats,
    curr_taxon_counts: &mut HashMap<TaxId, TaxonCounter>,
) {
    if call != 0 {
        stats.total_classified += 1;
        curr_taxon_counts
            .entry(call)
            .or_default()
            .increment_read_count();
    }
}

fn write_classification_output(
    call: TaxId,
    opts: &Options,
    taxonomy: &Taxonomy,
    dna: &Sequence,
    dna2: &Sequence,
    koss: &mut String,
    taxa: &[TaxId],
) {
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

    let ext_call = taxonomy.nodes()[&call].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            &taxonomy.name_data()[taxonomy.nodes()[&call].name_offset..]
        } else {
            "unclassified"
        };
        koss.push_str(name);
        koss.push_str(" (taxid ");
        koss.push_str(&ext_call.to_string());
        koss.push(')');
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
}

fn add_hitlist_string(oss: &mut String, taxa: &[TaxId], taxonomy: &Taxonomy) {
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in taxa.iter().skip(1) {
        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    oss.push_str(&format!("A:{} ", code_count));
                } else {
                    let ext_code = taxonomy.nodes[last_code as usize].external_id;
                    oss.push_str(&format!("{}:{} ", ext_code, code_count));
                }
            } else {
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
            let ext_code = taxonomy.nodes[last_code as usize].external_id;
            oss.push_str(&format!("{}:{}", ext_code, code_count));
        }
    } else {
        oss.push_str(if last_code == MATE_PAIR_BORDER_TAXON {
            "|:|"
        } else {
            "-:-"
        });
    }
}

fn initialize_outputs(opts: &Options, outputs: &mut OutputStreamData, format: SequenceFormat) {
    if !outputs.initialized {
        if !opts.classified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<_> = opts.classified_output_filename.split('#').collect();
                if fields.len() < 2 {
                    eprintln!(
                        "Paired filename format missing # character: {}",
                        opts.classified_output_filename
                    );
                    std::process::exit(1);
                } else if fields.len() > 2 {
                    eprintln!(
                        "Paired filename format has >1 # character: {}",
                        opts.classified_output_filename
                    );
                    std::process::exit(1);
                }
                outputs.classified_output1 = Some(BufWriter::new(
                    File::create(&format!("{}_1{}", fields[0], fields[1])).unwrap(),
                ));
                outputs.classified_output2 = Some(BufWriter::new(
                    File::create(&format!("{}_2{}", fields[0], fields[1])).unwrap(),
                ));
            } else {
                outputs.classified_output1 = Some(BufWriter::new(
                    File::create(&opts.classified_output_filename).unwrap(),
                ));
            }
            outputs.printing_sequences = true;
        }
        if !opts.unclassified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<_> = opts.unclassified_output_filename.split('#').collect();
                if fields.len() < 2 {
                    eprintln!(
                        "Paired filename format missing # character: {}",
                        opts.unclassified_output_filename
                    );
                    std::process::exit(1);
                } else if fields.len() > 2 {
                    eprintln!(
                        "Paired filename format has >1 # character: {}",
                        opts.unclassified_output_filename
                    );
                    std::process::exit(1);
                }
                outputs.unclassified_output1 = Some(BufWriter::new(
                    File::create(&format!("{}_1{}", fields[0], fields[1])).unwrap(),
                ));
                outputs.unclassified_output2 = Some(BufWriter::new(
                    File::create(&format!("{}_2{}", fields[0], fields[1])).unwrap(),
                ));
            } else {
                outputs.unclassified_output1 = Some(BufWriter::new(
                    File::create(&opts.unclassified_output_filename).unwrap(),
                ));
            }
            outputs.printing_sequences = true;
        }
        if !opts.kraken_output_filename.is_empty() {
            if opts.kraken_output_filename == "-" {
                outputs.kraken_output = None
            } else {
                outputs.kraken_output = Some(BufWriter::new(
                    File::create(&opts.kraken_output_filename).unwrap(),
                ));
            }
        }
        outputs.initialized = true;
    }
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: u8) {
    if dna.format != SequenceFormat::Fastq {
        return;
    }
    if dna.seq.len() != dna.quals.len() {
        eprintln!(
            "{}: Sequence length ({}) != Quality string length ({})",
            dna.id,
            dna.seq.len(),
            dna.quals.len()
        );
        std::process::exit(1);
    }
    for (base, qual) in dna.seq.as_bytes().iter_mut().zip(dna.quals.chars()) {
        if (qual as u8) - b'!' < minimum_quality_score {
            *base = b'x';
        }
    }
}
use clap::{App, Arg};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

fn parse_command_line(opts: &mut Options) {
    let matches = App::new("classify")
        .arg(
            Arg::with_name("index")
                .short("H")
                .long("index")
                .value_name("filename")
                .help("Kraken 2 index filename")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("taxonomy")
                .short("t")
                .long("taxonomy")
                .value_name("filename")
                .help("Kraken 2 taxonomy filename")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("options")
                .short("o")
                .long("options")
                .value_name("filename")
                .help("Kraken 2 options filename")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("quick")
                .short("q")
                .long("quick")
                .help("Quick mode"),
        )
        .arg(
            Arg::with_name("memory_mapping")
                .short("M")
                .long("memory-mapping")
                .help("Use memory mapping to access hash & taxonomy"),
        )
        .arg(
            Arg::with_name("confidence_threshold")
                .short("T")
                .long("confidence-threshold")
                .value_name("NUM")
                .help("Confidence score threshold (def. 0)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("threads")
                .short("p")
                .long("threads")
                .value_name("NUM")
                .help("Number of threads (def. 1)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("minimum_quality_score")
                .short("Q")
                .long("minimum-quality-score")
                .value_name("NUM")
                .help("Minimum quality score (FASTQ only, def. 0)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("paired")
                .short("P")
                .long("paired")
                .help("Process pairs of reads"),
        )
        .arg(
            Arg::with_name("single_file_pairs")
                .short("S")
                .long("single-file-pairs")
                .help("Process pairs with mates in same file"),
        )
        .arg(
            Arg::with_name("report")
                .short("R")
                .long("report")
                .value_name("filename")
                .help("Print report to filename")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mpa_style")
                .short("m")
                .long("mpa-style")
                .help("In comb. w/ -R, use mpa-style report"),
        )
        .arg(
            Arg::with_name("report_zero_counts")
                .short("z")
                .long("report-zero-counts")
                .help("In comb. w/ -R, report taxa w/ 0 count"),
        )
        .arg(
            Arg::with_name("scientific_name")
                .short("n")
                .long("scientific-name")
                .help("Print scientific name instead of taxid in Kraken output"),
        )
        .arg(
            Arg::with_name("minimum_hit_groups")
                .short("g")
                .long("minimum-hit-groups")
                .value_name("NUM")
                .help("Minimum number of hit groups needed for call")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("classified_output")
                .short("C")
                .long("classified-output")
                .value_name("filename")
                .help("Filename/format to have classified sequences")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("unclassified_output")
                .short("U")
                .long("unclassified-output")
                .value_name("filename")
                .help("Filename/format to have unclassified sequences")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("kraken_output")
                .short("O")
                .long("kraken-output")
                .value_name("filename")
                .help("Output file for normal Kraken output")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("kmer_data")
                .short("K")
                .long("kmer-data")
                .help("In comb. w/ -R, provide minimizer information in report"),
        )
        .get_matches();

    opts.index_filename = matches.value_of("index").unwrap().to_string();
    opts.taxonomy_filename = matches.value_of("taxonomy").unwrap().to_string();
    opts.options_filename = matches.value_of("options").unwrap().to_string();
    opts.quick_mode = matches.is_present("quick");
    opts.use_memory_mapping = matches.is_present("memory_mapping");

    if let Some(threshold) = matches.value_of("confidence_threshold") {
        opts.confidence_threshold = threshold.parse().unwrap();
        if opts.confidence_threshold < 0.0 || opts.confidence_threshold > 1.0 {
            eprintln!("confidence threshold must be in [0, 1]");
            std::process::exit(1);
        }
    }

    if let Some(threads) = matches.value_of("threads") {
        opts.num_threads = threads.parse().unwrap();
        if opts.num_threads < 1 {
            eprintln!("number of threads can't be less than 1");
            std::process::exit(1);
        }
    }

    if let Some(score) = matches.value_of("minimum_quality_score") {
        opts.minimum_quality_score = score.parse().unwrap();
    }

    opts.paired_end_processing = matches.is_present("paired");
    opts.single_file_pairs = matches.is_present("single_file_pairs");

    if let Some(filename) = matches.value_of("report") {
        opts.report_filename = filename.to_string();
    }

    opts.mpa_style_report = matches.is_present("mpa_style");
    opts.report_zero_counts = matches.is_present("report_zero_counts");
    opts.print_scientific_name = matches.is_present("scientific_name");

    if let Some(groups) = matches.value_of("minimum_hit_groups") {
        opts.minimum_hit_groups = groups.parse().unwrap();
    }

    if let Some(filename) = matches.value_of("classified_output") {
        opts.classified_output_filename = filename.to_string();
    }

    if let Some(filename) = matches.value_of("unclassified_output") {
        opts.unclassified_output_filename = filename.to_string();
    }

    if let Some(filename) = matches.value_of("kraken_output") {
        opts.kraken_output_filename = filename.to_string();
    }

    opts.report_kmer_data = matches.is_present("kmer_data");

    if opts.mpa_style_report && opts.report_filename.is_empty() {
        eprintln!("-m requires -R be used");
        usage(1);
    }
}

fn usage(exit_code: i32) {
    eprintln!("Usage: classify [options] <fasta/fastq file(s)>");
    eprintln!();
    eprintln!("Options: (*mandatory)");
    eprintln!("* -H, --index <filename>       Kraken 2 index filename");
    eprintln!("* -t, --taxonomy <filename>    Kraken 2 taxonomy filename");
    eprintln!("* -o, --options <filename>     Kraken 2 options filename");
    eprintln!("  -q, --quick                  Quick mode");
    eprintln!("  -M, --memory-mapping         Use memory mapping to access hash & taxonomy");
    eprintln!("  -T, --confidence-threshold <NUM>  Confidence score threshold (def. 0)");
    eprintln!("  -p, --threads <NUM>          Number of threads (def. 1)");
    eprintln!("  -Q, --minimum-quality-score <NUM>  Minimum quality score (FASTQ only, def. 0)");
    eprintln!("  -P, --paired                 Process pairs of reads");
    eprintln!("  -S, --single-file-pairs      Process pairs with mates in same file");
    eprintln!("  -R, --report <filename>      Print report to filename");
    eprintln!("  -m, --mpa-style              In comb. w/ -R, use mpa-style report");
    eprintln!("  -z, --report-zero-counts     In comb. w/ -R, report taxa w/ 0 count");
    eprintln!(
        "  -n, --scientific-name        Print scientific name instead of taxid in Kraken output"
    );
    eprintln!("  -g, --minimum-hit-groups <NUM>  Minimum number of hit groups needed for call");
    eprintln!("  -C, --classified-output <filename>  Filename/format to have classified sequences");
    eprintln!(
        "  -U, --unclassified-output <filename>  Filename/format to have unclassified sequences"
    );
    eprintln!("  -O, --kraken-output <filename>  Output file for normal Kraken output");
    eprintln!(
        "  -K, --kmer-data              In comb. w/ -R, provide minimizer information in report"
    );
    std::process::exit(exit_code);
}

fn main() {
    let mut opts = Options {
        quick_mode: false,
        confidence_threshold: 0.0,
        paired_end_processing: false,
        single_file_pairs: false,
        num_threads: 1,
        mpa_style_report: false,
        report_kmer_data: false,
        report_zero_counts: false,
        use_translated_search: false,
        print_scientific_name: false,
        minimum_quality_score: 0,
        minimum_hit_groups: 0,
        use_memory_mapping: false,
        index_filename: String::new(),
        taxonomy_filename: String::new(),
        options_filename: String::new(),
        report_filename: String::new(),
        classified_output_filename: todo!(),
        unclassified_output_filename: todo!(),
        kraken_output_filename: todo!(),
        match_input_order: todo!(),
    };

    let mut taxon_counters = HashMap::new();

    parse_command_line(&mut opts);

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    println!("Loading database information...");

    let idx_opts = load_index_options(&opts.options_filename);
    opts.use_translated_search = !idx_opts.dna_db;

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping);
    /// DEBUG: Added keybits/valuebits
    let hash_ptr = CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping, 111);

    println!("done.");

    let mut stats = ClassificationStats {
        total_sequences: 0,
        total_bases: 0,
        total_classified: 0,
    };

    let mut outputs = OutputStreamData {
        initialized: false,
        printing_sequences: false,
        kraken_output: None,
        classified_output1: None,
        classified_output2: None,
        unclassified_output1: None,
        unclassified_output2: Some(Box::new(std::io::stdout())),
    };

    let start_time = Instant::now();

    process_input_files(
        &opts,
        &idx_opts,
        &taxonomy,
        &hash_ptr,
        &mut stats,
        &mut outputs,
        &mut taxon_counters,
    );

    let end_time = Instant::now();

    report_stats(start_time, end_time, &stats);

    if !opts.report_filename.is_empty() {
        if opts.mpa_style_report {
            report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &taxon_counters,
            );
        } else {
            let total_unclassified = stats.total_sequences - stats.total_classified;
            report_kraken_style(
                &opts.report_filename,
                opts.report_zero_counts,
                opts.report_kmer_data,
                &taxonomy,
                &taxon_counters,
                stats.total_sequences,
                total_unclassified,
            );
        }
    }
}
