// SPDX-License-Identifier: MIT
// Author: Derrick Wood <dwood@cs.jhu.edu>

use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{IndexOptions, TaxId, TaxonCounters, TaxonCounts};
use crate::minimizer_scanner::MinimizerScanner;
use crate::seqreader::{BatchSequenceReader, Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::process;
use std::sync::{Arc, Mutex};
use std::time::Instant;

// Constants
const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
const MATE_PAIR_BORDER_TAXON: TaxidT = TAXID_MAX;
const READING_FRAME_BORDER_TAXON: TaxidT = TAXID_MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxidT = TAXID_MAX - 2;

// Type Aliases
type TaxidT = u32;
type TaxonCountsT = HashMap<TaxidT, u32>;
type TaxonCountersT<TaxonCounter> = HashMap<TaxidT, TaxonCounter>;

// Structs
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
    num_threads: i32,
    paired_end_processing: bool,
    single_file_pairs: bool,
    minimum_quality_score: i32,
    minimum_hit_groups: i32,
    use_memory_mapping: bool,
    match_input_order: bool,
}

#[derive(Default)]
struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

#[derive(Default)]
struct OutputStreamData {
    initialized: bool,
    printing_sequences: bool,
    classified_output1: Option<Arc<Mutex<Box<dyn Write>>>>,
    classified_output2: Option<Arc<Mutex<Box<dyn Write>>>>,
    unclassified_output1: Option<Arc<Mutex<Box<dyn Write>>>>,
    unclassified_output2: Option<Arc<Mutex<Box<dyn Write>>>>,
    kraken_output: Option<Arc<Mutex<Box<dyn Write>>>>,
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

// Function Prototypes
// fn parse_command_line(args: &[String], opts: &mut Options);
// fn usage(exit_code: i32);
// fn process_files(
//     filename1: Option<&str>,
//     filename2: Option<&str>,
//     hash: &KeyValueStore,
//     tax: &Taxonomy,
//     idx_opts: &IndexOptions,
//     opts: &Options,
//     stats: &mut ClassificationStats,
//     outputs: &mut OutputStreamData,
//     total_taxon_counters: &mut TaxonCountersT,
// );
// fn classify_sequence(
//     dna: &Sequence,
//     dna2: &Sequence,
//     koss: &mut String,
//     hash: &KeyValueStore,
//     tax: &Taxonomy,
//     idx_opts: &IndexOptions,
//     opts: &Options,
//     stats: &mut ClassificationStats,
//     scanner: &mut MinimizerScanner,
//     taxa: &mut Vec<taxid_t>,
//     hit_counts: &mut taxon_counts_t,
//     tx_frames: &mut Vec<String>,
//     my_taxon_counts: &mut TaxonCountersT,
// ) -> taxid_t;
// fn add_hitlist_string(oss: &mut String, taxa: &[taxid_t], taxonomy: &Taxonomy);
// fn resolve_tree(
//     hit_counts: &taxon_counts_t,
//     taxonomy: &Taxonomy,
//     total_minimizers: usize,
//     opts: &Options,
// ) -> taxid_t;
// fn report_stats(time1: Instant, time2: Instant, stats: &ClassificationStats);
// fn initialize_outputs(opts: &Options, outputs: &mut OutputStreamData, format: SequenceFormat);
// fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: i32);

fn main() {
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

    let mut taxon_counters: TaxonCountersT = HashMap::new(); // stats per taxon
    let args: Vec<String> = env::args().collect();
    parse_command_line(&args, &mut opts);

    // Multithreading is handled by the Rayon crate in Rust, which provides a parallel iterator interface
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads as usize)
        .build_global()
        .unwrap();

    eprintln!("Loading database information...");

    let idx_opts = IndexOptions::default();
    let idx_opt_fs = File::open(&opts.options_filename).expect("unable to open options file");
    let metadata = idx_opt_fs.metadata().expect("unable to get filesize");
    let opts_filesize = metadata.len();
    let mut buf_reader = BufReader::new(idx_opt_fs);
    let mut idx_opts_buffer = vec![0; opts_filesize as usize];
    buf_reader
        .read_exact(&mut idx_opts_buffer)
        .expect("unable to read options file");
    let idx_opts: IndexOptions =
        bincode::deserialize(&idx_opts_buffer).expect("unable to deserialize options");
    opts.use_translated_search = !idx_opts.dna_db;

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping);
    let hash_ptr = CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping);

    eprintln!(" done.");

    let mut stats = ClassificationStats::default();
    let mut outputs = OutputStreamData::default();

    let tv1 = Instant::now();
    if args.len() == 1 {
        if opts.paired_end_processing && !opts.single_file_pairs {
            eprintln!("paired end processing used with no files specified");
            process::exit(EX_USAGE);
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
        );
    } else {
        for i in 1..args.len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == args.len() {
                    eprintln!("paired end processing used with unpaired file");
                    process::exit(EX_USAGE);
                }
                process_files(
                    Some(&args[i]),
                    Some(&args[i + 1]),
                    &hash_ptr,
                    &taxonomy,
                    &idx_opts,
                    &opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                );
            } else {
                process_files(
                    Some(&args[i]),
                    None,
                    &hash_ptr,
                    &taxonomy,
                    &idx_opts,
                    &opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                );
            }
        }
    }
    let tv2 = Instant::now();

    drop(hash_ptr);

    report_stats(tv1, tv2, &stats);

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

fn report_stats(time1: Instant, time2: Instant, stats: &ClassificationStats) {
    let duration = time2.duration_since(time1);
    let seconds = duration.as_secs_f64();

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

fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: &KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    total_taxon_counters: &mut TaxonCountersT,
) {
    let fptr1: Box<dyn Read> = if let Some(filename) = filename1 {
        Box::new(File::open(filename).expect("unable to open file"))
    } else {
        Box::new(io::stdin())
    };

    let fptr2: Option<Box<dyn Read>> = if let Some(filename) = filename2 {
        Some(Box::new(File::open(filename).expect("unable to open file")))
    } else {
        None
    };

    // The priority queue for output is designed to ensure fragment data
    // is output in the same order it was input
    let output_queue: Arc<Mutex<Vec<OutputData>>> = Arc::new(Mutex::new(Vec::new()));
    let next_input_block_id = Arc::new(Mutex::new(0));
    let next_output_block_id = Arc::new(Mutex::new(0));

    rayon::scope(|s| {
        s.spawn(|_| {
            let mut scanner = MinimizerScanner::new(
                idx_opts.k,
                idx_opts.l,
                idx_opts.spaced_seed_mask,
                idx_opts.dna_db,
                idx_opts.toggle_mask,
                idx_opts.revcom_version,
            );
            let mut taxa: Vec<TaxidT> = Vec::new();
            let mut hit_counts: TaxonCountsT = HashMap::new();
            let mut kraken_oss = String::new();
            let mut c1_oss = String::new();
            let mut c2_oss = String::new();
            let mut u1_oss = String::new();
            let mut u2_oss = String::new();
            let mut thread_stats = ClassificationStats::default();
            let mut translated_frames: Vec<String> = vec![String::new(); 6];
            let mut reader1 = BatchSequenceReader::new();
            let mut reader2 = BatchSequenceReader::new();
            use crate::sequence::Sequence;

            let mut seq1 = Sequence::default();
            let mut seq2 = Sequence::default();
            let mut thread_taxon_counters: TaxonCountersT = HashMap::new();

            loop {
                thread_stats.total_sequences = 0;
                thread_stats.total_bases = 0;
                thread_stats.total_classified = 0;

                let mut ok_read = false;

                {
                    let mut next_input_block_id = next_input_block_id.lock().unwrap();
                    // Input processing block
                    if !opts.paired_end_processing {
                        // Unpaired data? Just read in a sized block
                        ok_read = reader1.load_block(&fptr1, 3 * 1024 * 1024);
                    } else if !opts.single_file_pairs {
                        // Paired data in 2 files? Read a line-counted batch from each file.
                        ok_read = reader1.load_batch(&fptr1, NUM_FRAGMENTS_PER_THREAD);
                        if ok_read && opts.paired_end_processing {
                            ok_read = reader2.load_batch(&fptr2.unwrap(), NUM_FRAGMENTS_PER_THREAD);
                        }
                    } else {
                        let frags = NUM_FRAGMENTS_PER_THREAD * 2;
                        // Ensure frag count is even - just in case above line is changed
                        ok_read = reader1.load_batch(&fptr1, frags);
                    }
                    *next_input_block_id += 1;
                }

                if !ok_read {
                    break;
                }

                // Reset all dynamically-growing things
                kraken_oss.clear();
                c1_oss.clear();
                c2_oss.clear();
                u1_oss.clear();
                u2_oss.clear();
                thread_taxon_counters.clear();

                while let Some(valid_fragment) = reader1.next_sequence(&mut seq1) {
                    if opts.paired_end_processing && valid_fragment {
                        if opts.single_file_pairs {
                            reader1.next_sequence(&mut seq2);
                        } else {
                            reader2.next_sequence(&mut seq2);
                        }
                    }
                    if !valid_fragment {
                        break;
                    }
                    thread_stats.total_sequences += 1;
                    if opts.minimum_quality_score > 0 {
                        mask_low_quality_bases(&mut seq1, opts.minimum_quality_score);
                        if opts.paired_end_processing {
                            mask_low_quality_bases(&mut seq2, opts.minimum_quality_score);
                        }
                    }
                    let call = classify_sequence(
                        &seq1,
                        &seq2,
                        &mut kraken_oss,
                        &hash,
                        &tax,
                        &idx_opts,
                        &opts,
                        &mut thread_stats,
                        &mut scanner,
                        &mut taxa,
                        &mut hit_counts,
                        &mut translated_frames,
                        &mut thread_taxon_counters,
                    );
                    if call != 0 {
                        let buffer = format!(" kraken:taxid|{}", tax.nodes()[call].external_id);
                        seq1.header.push_str(&buffer);
                        c1_oss.push_str(&seq1.to_string());
                        if opts.paired_end_processing {
                            seq2.header.push_str(&buffer);
                            c2_oss.push_str(&seq2.to_string());
                        }
                    } else {
                        u1_oss.push_str(&seq1.to_string());
                        if opts.paired_end_processing {
                            u2_oss.push_str(&seq2.to_string());
                        }
                    }
                    thread_stats.total_bases += seq1.seq.len() as u64;
                    if opts.paired_end_processing {
                        thread_stats.total_bases += seq2.seq.len() as u64;
                    }
                }

                stats.total_sequences += thread_stats.total_sequences;
                stats.total_bases += thread_stats.total_bases;
                stats.total_classified += thread_stats.total_classified;

                if io::stderr().is_terminal() {
                    eprintln!(
                        "\rProcessed {} sequences ({} bp) ...",
                        stats.total_sequences, stats.total_bases
                    );
                }

                if !outputs.initialized {
                    initialize_outputs(&opts, outputs, reader1.file_format());
                }

                let out_data = OutputData {
                    block_id: *next_input_block_id.lock().unwrap(),
                    kraken_str: kraken_oss.clone(),
                    classified_out1_str: c1_oss.clone(),
                    classified_out2_str: c2_oss.clone(),
                    unclassified_out1_str: u1_oss.clone(),
                    unclassified_out2_str: u2_oss.clone(),
                };

                output_queue.lock().unwrap().push(out_data);

                for (key, value) in thread_taxon_counters.iter() {
                    *total_taxon_counters.entry(*key).or_insert(0) += *value;
                }

                loop {
                    let mut output_queue = output_queue.lock().unwrap();
                    if output_queue.is_empty()
                        || output_queue[0].block_id != *next_output_block_id.lock().unwrap()
                    {
                        break;
                    }

                    let out_data = output_queue.remove(0);
                    if let Some(ref kraken_output) = outputs.kraken_output {
                        let mut kraken_output = kraken_output.lock().unwrap();
                        kraken_output
                            .write_all(out_data.kraken_str.as_bytes())
                            .unwrap();
                    }
                    if let Some(ref classified_output1) = outputs.classified_output1 {
                        let mut classified_output1 = classified_output1.lock().unwrap();
                        classified_output1
                            .write_all(out_data.classified_out1_str.as_bytes())
                            .unwrap();
                    }
                    if let Some(ref classified_output2) = outputs.classified_output2 {
                        let mut classified_output2 = classified_output2.lock().unwrap();
                        classified_output2
                            .write_all(out_data.classified_out2_str.as_bytes())
                            .unwrap();
                    }
                    if let Some(ref unclassified_output1) = outputs.unclassified_output1 {
                        let mut unclassified_output1 = unclassified_output1.lock().unwrap();
                        unclassified_output1
                            .write_all(out_data.unclassified_out1_str.as_bytes())
                            .unwrap();
                    }
                    if let Some(ref unclassified_output2) = outputs.unclassified_output2 {
                        let mut unclassified_output2 = unclassified_output2.lock().unwrap();
                        unclassified_output2
                            .write_all(out_data.unclassified_out2_str.as_bytes())
                            .unwrap();
                    }

                    *next_output_block_id.lock().unwrap() += 1;
                }
            }
        });
    });

    if let Some(fptr1) = fptr1 {
        drop(fptr1);
    }
    if let Some(fptr2) = fptr2 {
        drop(fptr2);
    }
    if let Some(ref kraken_output) = outputs.kraken_output {
        let mut kraken_output = kraken_output.lock().unwrap();
        kraken_output.flush().unwrap();
    }
    if let Some(ref classified_output1) = outputs.classified_output1 {
        let mut classified_output1 = classified_output1.lock().unwrap();
        classified_output1.flush().unwrap();
    }
    if let Some(ref classified_output2) = outputs.classified_output2 {
        let mut classified_output2 = classified_output2.lock().unwrap();
        classified_output2.flush().unwrap();
    }
    if let Some(ref unclassified_output1) = outputs.unclassified_output1 {
        let mut unclassified_output1 = unclassified_output1.lock().unwrap();
        unclassified_output1.flush().unwrap();
    }
    if let Some(ref unclassified_output2) = outputs.unclassified_output2 {
        let mut unclassified_output2 = unclassified_output2.lock().unwrap();
        unclassified_output2.flush().unwrap();
    }
}

fn classify_sequence(
    dna: &Sequence,
    dna2: &Sequence,
    koss: &mut String,
    hash: &KeyValueStore,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<TaxidT>,
    hit_counts: &mut TaxonCountsT,
    tx_frames: &mut Vec<String>,
    curr_taxon_counts: &mut TaxonCountersT,
) -> TaxidT {
    let mut minimizer_ptr: Option<&u64>;
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
            translate_to_all_frames(if mate_num == 0 { &dna.seq } else { &dna2.seq }, tx_frames);
        }

        for frame_idx in 0..frame_ct {
            if opts.use_translated_search {
                scanner.load_sequence(&tx_frames[frame_idx]);
            } else {
                scanner.load_sequence(if mate_num == 0 { &dna.seq } else { &dna2.seq });
            }
            let mut last_minimizer = u64::MAX;
            let mut last_taxon = TAXID_MAX;
            while {
                minimizer_ptr = scanner.next_minimizer();
                minimizer_ptr.is_some()
            } {
                let taxon;
                if scanner.is_ambiguous() {
                    taxon = AMBIGUOUS_SPAN_TAXON;
                } else {
                    let minimizer = *minimizer_ptr.unwrap();
                    if minimizer != last_minimizer {
                        let mut skip_lookup = false;
                        if let Some(min_acceptable_hash) = idx_opts.minimum_acceptable_hash_value {
                            if murmur_hash3(minimizer) < min_acceptable_hash {
                                skip_lookup = true;
                            }
                        }
                        taxon = if skip_lookup { 0 } else { hash.get(minimizer) };
                        last_taxon = taxon;
                        last_minimizer = minimizer;
                        if taxon != 0 {
                            minimizer_hit_groups += 1;
                            curr_taxon_counts
                                .entry(taxon)
                                .or_insert_with(|| TaxonCounter::new())
                                .add_kmer(scanner.last_minimizer());
                        }
                    } else {
                        taxon = last_taxon;
                    }
                    if taxon != 0 {
                        if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                            call = taxon;
                            break;
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

    if call == 0 {
        let total_kmers = taxa.len()
            - if opts.paired_end_processing { 1 } else { 0 }
            - if opts.use_translated_search {
                if opts.paired_end_processing {
                    4
                } else {
                    2
                }
            } else {
                0
            };
        call = resolve_tree(hit_counts, taxonomy, total_kmers, opts);
        if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
            call = 0;
        }
    }

    if call != 0 {
        stats.total_classified += 1;
        curr_taxon_counts
            .entry(call)
            .or_insert_with(|| TaxonCounter::new())
            .increment_read_count();
    }

    koss.push_str(if call != 0 { "C\t" } else { "U\t" });
    koss.push_str(&format!(
        "{}\t",
        if !opts.paired_end_processing {
            &dna.id
        } else {
            &trim_pair_info(&dna.id)
        }
    ));

    let ext_call = taxonomy.nodes()[call].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            Some(taxonomy.name_data() + taxonomy.nodes()[call].name_offset)
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

    koss.push_str("\t");
    koss.push_str(&format!(
        "{}\t",
        if !opts.paired_end_processing {
            dna.seq.len().to_string()
        } else {
            format!("{}|{}", dna.seq.len(), dna2.seq.len())
        }
    ));

    if opts.quick_mode {
        koss.push_str(&format!("{}:Q", ext_call));
    } else {
        if taxa.is_empty() {
            koss.push_str("0:0");
        } else {
            add_hitlist_string(koss, &taxa, taxonomy);
        }
    }

    koss.push('\n');

    call
}

fn add_hitlist_string(oss: &mut String, taxa: &[TaxidT], taxonomy: &Taxonomy) {
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in &taxa[1..] {
        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    oss.push_str(&format!("A:{} ", code_count));
                } else {
                    let ext_code = taxonomy.nodes()[last_code].external_id;
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
            let ext_code = taxonomy.nodes()[last_code].external_id;
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
                let fields: Vec<&str> = opts.classified_output_filename.splitn(3, '#').collect();
                if fields.len() < 2 {
                    eprintln!(
                        "Paired filename format missing # character: {}",
                        opts.classified_output_filename
                    );
                    process::exit(EX_DATAERR);
                } else if fields.len() > 2 {
                    eprintln!(
                        "Paired filename format has >1 # character: {}",
                        opts.classified_output_filename
                    );
                    process::exit(EX_DATAERR);
                }
                outputs.classified_output1 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(format!("{}_1{}", fields[0], fields[1]))
                        .unwrap(),
                ))));
                outputs.classified_output2 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(format!("{}_2{}", fields[0], fields[1]))
                        .unwrap(),
                ))));
            } else {
                outputs.classified_output1 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(&opts.classified_output_filename)
                        .unwrap(),
                ))));
            }
            outputs.printing_sequences = true;
        }
        if !opts.unclassified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<&str> = opts.unclassified_output_filename.splitn(3, '#').collect();
                if fields.len() < 2 {
                    eprintln!(
                        "Paired filename format missing # character: {}",
                        opts.unclassified_output_filename
                    );
                    process::exit(EX_DATAERR);
                } else if fields.len() > 2 {
                    eprintln!(
                        "Paired filename format has >1 # character: {}",
                        opts.unclassified_output_filename
                    );
                    process::exit(EX_DATAERR);
                }
                outputs.unclassified_output1 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(format!("{}_1{}", fields[0], fields[1]))
                        .unwrap(),
                ))));
                outputs.unclassified_output2 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(format!("{}_2{}", fields[0], fields[1]))
                        .unwrap(),
                ))));
            } else {
                outputs.unclassified_output1 = Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(&opts.unclassified_output_filename)
                        .unwrap(),
                ))));
            }
            outputs.printing_sequences = true;
        }
        if !opts.kraken_output_filename.is_empty() {
            outputs.kraken_output = if opts.kraken_output_filename == "-" {
                None
            } else {
                Some(Arc::new(Mutex::new(Box::new(
                    OpenOptions::new()
                        .create(true)
                        .write(true)
                        .open(&opts.kraken_output_filename)
                        .unwrap(),
                ))))
            };
        }
        outputs.initialized = true;
    }
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: i32) {
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
        process::exit(EX_DATAERR);
    }
    for (i, &qual) in dna.quals.iter().enumerate() {
        if (qual as i32 - '!' as i32) < minimum_quality_score {
            dna.seq[i] = 'x';
        }
    }
}

fn parse_command_line(args: &[String], opts: &mut Options) {
    let mut optind = 1;
    while optind < args.len() {
        match args[optind].as_str() {
            "-h" | "?" => usage(0),
            "-H" => opts.index_filename = args[optind + 1].clone(),
            "-t" => opts.taxonomy_filename = args[optind + 1].clone(),
            "-T" => {
                opts.confidence_threshold = args[optind + 1]
                    .parse()
                    .expect("confidence threshold must be a number");
                if opts.confidence_threshold < 0.0 || opts.confidence_threshold > 1.0 {
                    eprintln!("confidence threshold must be in [0, 1]");
                    process::exit(EX_USAGE);
                }
            }
            "-o" => opts.options_filename = args[optind + 1].clone(),
            "-q" => opts.quick_mode = true,
            "-p" => {
                opts.num_threads = args[optind + 1]
                    .parse()
                    .expect("number of threads must be an integer");
                if opts.num_threads < 1 {
                    eprintln!("number of threads can't be less than 1");
                    process::exit(EX_USAGE);
                }
            }
            "-g" => {
                opts.minimum_hit_groups = args[optind + 1]
                    .parse()
                    .expect("minimum hit groups must be an integer")
            }
            "-P" => opts.paired_end_processing = true,
            "-S" => {
                opts.paired_end_processing = true;
                opts.single_file_pairs = true;
            }
            "-m" => opts.mpa_style_report = true,
            "-K" => opts.report_kmer_data = true,
            "-R" => opts.report_filename = args[optind + 1].clone(),
            "-z" => opts.report_zero_counts = true,
            "-C" => opts.classified_output_filename = args[optind + 1].clone(),
            "-U" => opts.unclassified_output_filename = args[optind + 1].clone(),
            "-O" => opts.kraken_output_filename = args[optind + 1].clone(),
            "-n" => opts.print_scientific_name = true,
            "-Q" => {
                opts.minimum_quality_score = args[optind + 1]
                    .parse()
                    .expect("minimum quality score must be an integer")
            }
            "-M" => opts.use_memory_mapping = true,
            _ => {
                eprintln!("Unknown option: {}", args[optind]);
                usage(EX_USAGE);
            }
        }
        optind += 2;
    }

    if opts.index_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        eprintln!("mandatory filename missing");
        usage(EX_USAGE);
    }

    if opts.mpa_style_report && opts.report_filename.is_empty() {
        eprintln!("-m requires -R be used");
        usage(EX_USAGE);
    }
}

fn usage(exit_code: i32) {
    eprintln!(
        "Usage: classify [options] <fasta/fastq file(s)>
Options: (*mandatory)
* -H filename      Kraken 2 index filename
* -t filename      Kraken 2 taxonomy filename
* -o filename      Kraken 2 options filename
  -q               Quick mode
  -M               Use memory mapping to access hash & taxonomy
  -T NUM           Confidence score threshold (def. 0)
  -p NUM           Number of threads (def. 1)
  -Q NUM           Minimum quality score (FASTQ only, def. 0)
  -P               Process pairs of reads
  -S               Process pairs with mates in same file
  -R filename      Print report to filename
  -m               In comb. w/ -R, use mpa-style report
  -z               In comb. w/ -R, report taxa w/ 0 count
  -n               Print scientific name instead of taxid in Kraken output
  -g NUM           Minimum number of hit groups needed for call
  -C filename      Filename/format to have classified sequences
  -U filename      Filename/format to have unclassified sequences
  -O filename      Output file for normal Kraken output
  -K               In comb. w/ -R, provide minimizer information in report"
    );
    process::exit(exit_code);
}
