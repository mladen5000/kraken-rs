use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{IndexOptions, TaxId, TAXID_MAX};
use crate::kv_store::{self, KeyValueStore};
use crate::mmscanner::MinimizerScanner;
use crate::reports::{self, TaxonCounters, TaxonCounts};
use crate::seqreader::{Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy;

use clap::{command, Parser};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::fs::{metadata, File};
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::{Arc, Mutex};
use std::time::Instant;

pub const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
pub const MATE_PAIR_BORDER_TAXON: TaxId = TAXID_MAX;
pub const READING_FRAME_BORDER_TAXON: TaxId = TAXID_MAX - 1;
pub const AMBIGUOUS_SPAN_TAXON: TaxId = TAXID_MAX - 2;

/// Application options specified via command-line arguments.
#[command(version, about, long_about = None)]
#[derive(Debug, Parser)]
struct Options {
    /// Kraken 2 index filename
    #[clap(short = 'H', long)]
    index_filename: String,

    /// Kraken 2 taxonomy filename
    #[clap(short = 't', long)]
    taxonomy_filename: String,

    /// Kraken 2 options filename
    #[clap(short = 'o', long)]
    options_filename: String,

    /// Filename/format to have classified sequences
    #[clap(short = 'C', long)]
    classified_output_filename: String,

    /// Filename/format to have unclassified sequences
    #[clap(short = 'U', long)]
    unclassified_output_filename: String,

    /// Output file for normal Kraken output
    #[clap(short = 'O', long)]
    kraken_output_filename: String,

    /// Quick mode
    #[clap(short = 'q', long, action)]
    quick_mode: bool,

    /// Use memory mapping to access hash & taxonomy
    #[clap(short = 'M', long, action)]
    use_memory_mapping: bool,

    /// Confidence score threshold
    #[clap(short = 'T', long, default_value_t = 0.0)]
    confidence_threshold: f32,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = 1)]
    num_threads: usize,

    /// Minimum quality score (FASTQ only)
    #[clap(short = 'Q', long, default_value_t = 0)]
    minimum_quality_score: u32,

    /// Process pairs of reads
    #[clap(short = 'P', long, action)]
    paired_end_processing: bool,

    /// Process pairs with mates in same file
    #[clap(short = 'S', long, action)]
    single_file_pairs: bool,

    /// Print report to filename
    #[clap(short = 'R', long)]
    report_filename: String,

    /// In combination with -R, use mpa-style report
    #[clap(short = 'm', long, action)]
    mpa_style_report: bool,

    /// In combination with -R, report taxa with 0 count
    #[clap(short = 'z', long, action)]
    report_zero_counts: bool,

    /// Print scientific name instead of taxid in Kraken output
    #[clap(short = 'n', long, action)]
    print_scientific_name: bool,

    /// Minimum number of hit groups needed for call
    #[clap(short = 'g', long, default_value_t = 0)]
    minimum_hit_groups: u32,

    /// In combination with -R, provide minimizer information in report
    #[clap(short = 'K', long, action)]
    report_kmer_data: bool,

    /// Use translated search
    #[clap(short = 'X', long, action)]
    use_translated_search: bool,

    #[clap(short = 'I', long, action)]
    match_input_order: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            index_filename: String::new(),
            taxonomy_filename: String::new(),
            options_filename: String::new(),
            report_filename: String::new(),
            classified_output_filename: String::new(),
            unclassified_output_filename: String::new(),
            kraken_output_filename: String::new(),
            mpa_style_report: false,
            report_kmer_data: false,
            quick_mode: false,
            report_zero_counts: false,
            print_scientific_name: false,
            confidence_threshold: 0.0,
            num_threads: 1,
            paired_end_processing: false,
            single_file_pairs: false,
            minimum_quality_score: 0,
            minimum_hit_groups: 0,
            use_memory_mapping: false,
            use_translated_search: false,
            match_input_order: false,
        }
    }
}

#[derive(Default)]
pub struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

#[derive(Default)]
pub struct OutputStreamData {
    classified_output1: Option<File>,
    classified_output2: Option<File>,
    unclassified_output1: Option<File>,
    unclassified_output2: Option<File>,
    kraken_output: Option<File>,
    initialized: bool,
    printing_sequences: bool,
}
#[derive(Default)]
pub struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}

fn main() {
    let opts = Options::default();
    let taxon_counters = TaxonCounters::new();

    parse_command_line(&mut opts);

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()
        .unwrap();

    eprintln!("Loading database information ...");

    let mut idx_opts = IndexOptions::new();

    let mut idx_opts_fs =
        File::open(&opts.options_filename).expect("Unable to get index options file");
    let opts_filesize = metadata(&opts.options_filename)
        .expect("Unable to get metadata")
        .len();
    let mut buffer = vec![0; opts_filesize as usize];
    idx_opts_fs
        .read_exact(&mut buffer)
        .expect("Unable to read file");
    idx_opts = unsafe { std::ptr::read(buffer.as_ptr() as *const _) };

    opts.use_translated_search = !idx_opts.dna_db;

    let mut taxonomy = Taxonomy::new(opts.taxonomy_filename.as_str(), opts.use_memory_mapping);
    let mut hash_ptr =
        CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping, &idx_opts);

    println!(" done.");

    let mut stats = ClassificationStats::default();

    let mut outputs = OutputStreamData::default();

    let start_time = Instant::now();
    if std::env::args().len() == 1 {
        if opts.paired_end_processing && !opts.single_file_pairs {
            panic!("paired end processing used with no files specified");
        }
        process_files(None, None, opts);
    } else {
        for i in 1..std::env::args().len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == std::env::args().len() {
                    panic!("paired end processing used with unpaired file");
                }
                let hash_ptr = Arc::new(Mutex::new(hash_ptr.unwrap()));
                process_files(
                    Some(&std::env::args().nth(i).unwrap()),
                    Some(&std::env::args().nth(i + 1).unwrap()),
                    opts,
                );
            } else {
                process_files(Some(&std::env::args().nth(i).unwrap()), None, opts);
            }
        }
    }

    let end_time = Instant::now();

    report_stats(start_time, end_time, &mut stats);

    let file = std::fs::File::create(&opts.report_filename).expect("Unable to create file");
    if !opts.report_filename.unwrap().is_empty() {
        if opts.mpa_style_report {
            reports::report_mpa_style(
                &mut file,
                opts.report_zero_counts,
                &mut taxonomy.unwrap(),
                &mut taxon_counters,
            );
        } else {
            let total_unclassified = stats.total_sequences - stats.total_classified;
            reports::report_kraken_style(
                &mut file,
                opts.report_zero_counts,
                opts.report_kmer_data,
                &mut taxonomy.unwrap(),
                &mut taxon_counters,
                stats.total_sequences,
                total_unclassified,
            );
        }
    }
}

fn report_stats(start_time: Instant, end_time: Instant, stats: &mut ClassificationStats) {
    let duration = end_time.duration_since(start_time);
    let seconds = duration.as_secs() as f64 + duration.subsec_nanos() as f64 * 1e-9;

    let total_unclassified = stats.total_sequences - stats.total_classified;

    println!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1.0e6,
        seconds,
        stats.total_sequences as f64 / 1.0e3 / (seconds / 60.0),
        stats.total_bases as f64 / 1.0e6 / (seconds / 60.0)
    );
    println!(
        "  {} sequences classified ({:.2}%)",
        stats.total_classified,
        stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
    );
    println!(
        "  {} sequences unclassified ({:.2}%)",
        total_unclassified,
        total_unclassified as f64 * 100.0 / stats.total_sequences as f64
    );
}

// fn process_files(
//     filename1: Option<&str>,
//     filename2: Option<&str>,
//     hash: Arc<Mutex<CompactHashTable>>,
//     tax: Arc<Mutex<Taxonomy>>,
//     idx_opts: Arc<Mutex<IndexOptions>>,
//     opts: Arc<Mutex<Options>>,
//     stats: Arc<Mutex<ClassificationStats>>,
//     outputs: Arc<Mutex<OutputStreamData>>,
//     total_taxon_counters: Arc<Mutex<TaxonCounters>>,
// ) {
//     let fptr1 = match filename1 {
//         Some(name) => Arc::new(Mutex::new(
//             Box::new(BufReader::new(File::open(name).unwrap())) as Box<dyn Read>,
//         )),
//         None => Arc::new(Mutex::new(
//             Box::new(BufReader::new(std::io::stdin())) as Box<dyn Read>
//         )),
//     };

//     let fptr2: Box<dyn Read> = if opts.paired_end_processing && !opts.single_file_pairs {
//         match filename2 {
//             Some(filename) => Box::new(BufReader::new(File::open(filename)?)),
//             None => {
//                 return Err(io::Error::new(
//                     io::ErrorKind::Other,
//                     "filename2 is required",
//                 ))
//             }
//         }
//     } else {
//         fptr1
//     };

//     let output_queue: Arc<Mutex<BinaryHeap<OutputData>>> = Arc::new(Mutex::new(BinaryHeap::new()));

//     let threads: Vec<_> = (0..num_cpus::get())
//         .map(|_| {
//             let queue_clone = Arc::clone(&output_queue);
//             let fptr1_clone = Arc::clone(&fptr1);
//             let fptr2_clone = fptr2.as_ref().map(Arc::clone);

//             thread::spawn(move || {
//                 // Thread-local variables go here

//                 loop {
//                     // Read and process a block of data from fptr1 and possibly fptr2
//                     let mut thread_stats = ClassificationStats {
//                         total_sequences: 0,
//                         total_bases: 0,
//                         total_classified: 0,
//                     };
//                     let mut kraken_oss = String::new();
//                     let mut c1_oss = String::new();
//                     let mut c2_oss = String::new();
//                     let mut u1_oss = String::new();
//                     let mut u2_oss = String::new();
//                     let mut thread_taxon_counters = HashMap::new();

//                     let reader1 = BufReader::new(fptr1.unwrap());
//                     let reader2 = fptr2.map(|fptr| BufReader::new(fptr.unwrap()));
//                     if !opts.paired_end_processing {
//                         // Unpaired data? Just read in a sized block
//                         ok_read = reader1.load_block(fptr1, 3 * 1024 * 1024);
//                     } else if !opts.single_file_pairs {
//                         // Paired data in 2 files? Read a line-counted batch from each file.
//                         ok_read = reader1.load_batch(fptr1, NUM_FRAGMENTS_PER_THREAD);
//                         if ok_read && opts.paired_end_processing {
//                             ok_read = reader2.load_batch(fptr2, NUM_FRAGMENTS_PER_THREAD);
//                         }
//                     } else {
//                         let mut frags = NUM_FRAGMENTS_PER_THREAD * 2;
//                         // Ensure frag count is even - just in case above line is changed
//                         if frags % 2 == 1 {
//                             frags += 1;
//                         }
//                         ok_read = reader1.load_batch(fptr1, frags);
//                     }
//                     let mut next_input_block_id = 0;

//                     // block_id = next_input_block_id.fetch_add(1, Ord::SeqCst);

//                     if !ok_read {
//                         break;
//                     }

//                     // Reset all dynamically-growing things
//                     kraken_oss.clear();
//                     c1_oss.clear();
//                     c2_oss.clear();
//                     u1_oss.clear();
//                     u2_oss.clear();
//                     thread_taxon_counters.clear();

//                     while let Some(mut seq1) = reader1.next_sequence() {
//                         let mut seq2 = if opts.paired_end_processing {
//                             if opts.single_file_pairs {
//                                 reader1.next_sequence()
//                             } else {
//                                 reader2.next_sequence()
//                             }
//                         } else {
//                             None
//                         };

//                         if seq2.is_none() {
//                             break;
//                         }

//                         thread_stats.total_sequences += 1;

//                         if opts.minimum_quality_score > 0 {
//                             mask_low_quality_bases(&mut seq1, opts.minimum_quality_score);
//                             if let Some(seq2) = &mut seq2 {
//                                 mask_low_quality_bases(seq2, opts.minimum_quality_score);
//                             }
//                         }

//                         let mut translated_frames = Vec::new();
//                         let mut hit_counts = HashMap::new(); // Declare the missing hit_counts variable
//                         let mut taxmer;
//                         let call = classify_sequence(
//                             &mut seq1,
//                             &mut seq2,
//                             &mut kraken_oss,
//                             &mut hash,
//                             &mut tax,
//                             &idx_opts,
//                             &opts,
//                             &mut thread_stats,
//                             &mut taxmer,
//                             &mut tax,
//                             &mut hit_counts, // Pass the hit_counts variable
//                             &mut translated_frames,
//                             &mut thread_taxon_counters,
//                         );

//                         if let Some(call) = call {
//                             let buffer = format!(" kraken:taxid|{}", tax.nodes()[call].external_id);
//                             seq1.header.push_str(&buffer);
//                             write!(&mut c1_oss, "{}", format!("{}", seq1)).unwrap();
//                             if let Some(mut seq2) = seq2 {
//                                 seq2.header.push_str(&buffer);
//                                 write!(&mut c2_oss, "{}", format!("{}", seq2)).unwrap();
//                             }
//                         } else {
//                             write!(&mut u1_oss, "{}", format!("{}", seq1)).unwrap();
//                             if let Some(seq2) = seq2 {
//                                 write!(&mut u2_oss, "{}", format!("{}", seq2)).unwrap();
//                             }
//                         }
//                         thread_stats.total_bases += seq1.seq.len();
//                         if let Some(seq2) = seq2 {
//                             thread_stats.total_bases += seq2.seq.len();
//                         }
//                     }

//                     stats
//                         .total_sequences
//                         .fetch_add(thread_stats.total_sequences, Ord::SeqCst);
//                     stats
//                         .total_bases
//                         .fetch_add(thread_stats.total_bases, Ord::SeqCst);
//                     stats
//                         .total_classified
//                         .fetch_add(thread_stats.total_classified, Ord::SeqCst);

//                     if atty::is(atty::Stream::Stderr) {
//                         eprint!(
//                             "\rProcessed {} sequences ({} bp) ...",
//                             stats.total_sequences.load(Ord::SeqCst),
//                             stats.total_bases.load(Ord::SeqCst)
//                         );
//                     }

//                     if !outputs.initialized {
//                         initialize_outputs(&opts, &mut outputs);
//                     }

//                     let block_id = 0; // Declare the block_id variable

//                     let mut out_data = OutputData {
//                         block_id,
//                         kraken_str: kraken_oss.clone(),
//                         classified_out1_str: c1_oss.clone(),
//                         classified_out2_str: c2_oss.clone(),
//                         unclassified_out1_str: u1_oss.clone(),
//                         unclassified_out2_str: u2_oss.clone(),
//                     };

//                     {
//                         let mut output_queue_guard = output_queue.lock().unwrap();
//                         output_queue_guard.push(out_data);
//                     }

//                     {
//                         let mut total_taxon_counters_guard = total_taxon_counters.lock().unwrap();
//                         for (key, value) in thread_taxon_counters.drain() {
//                             *total_taxon_counters_guard.entry(key).or_insert(0) += value;
//                         }
//                     }

//                     let mut output_loop = true;

//                     while output_loop {
//                         let mut next_output_block_id_guard = next_input_block_id.lock().unwrap();
//                         let mut output_queue_guard = output_queue.lock().unwrap();
//                         output_loop = !output_queue_guard.is_empty();
//                         if output_loop {
//                             out_data = output_queue_guard.peek().unwrap().clone();
//                             if out_data.block_id == *next_output_block_id_guard {
//                                 output_queue_guard.pop();
//                                 *next_output_block_id_guard += 1;
//                             } else {
//                                 output_loop = false;
//                             }
//                         }
//                         drop(next_output_block_id_guard);
//                         drop(output_queue_guard);

//                         if !output_loop {
//                             break;
//                         }

//                         if let Some(kraken_output) = &mut outputs.kraken_output {
//                             write!(kraken_output, "{}", out_data.kraken_str).unwrap();
//                         }
//                         if let Some(classified_output1) = &mut outputs.classified_output1 {
//                             write!(classified_output1, "{}", out_data.classified_out1_str).unwrap();
//                         }
//                         if let Some(classified_output2) = &mut outputs.classified_output2 {
//                             write!(classified_output2, "{}", out_data.classified_out2_str).unwrap();
//                         }
//                         if let Some(unclassified_output1) = &mut outputs.unclassified_output1 {
//                             write!(unclassified_output1, "{}", out_data.unclassified_out1_str)
//                                 .unwrap();
//                         }
//                     }
//                 } // Add closing parenthesis for the closure
//             })
//         })
//         .collect();
// }

// Flush the output streams

fn resolve_tree(
    hit_counts: &mut TaxonCounts,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (opts.confidence_threshold * total_minimizers as f32).ceil() as u32;

    for (&taxon, _) in hit_counts.iter() {
        let mut score = 0;

        for (&taxon2, &count) in hit_counts.iter() {
            if taxonomy.is_a_ancestor_of_b(taxon2 as usize, taxon as usize) {
                score += count;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score {
            max_taxon = taxonomy.lowest_common_ancestor(max_taxon as u64, taxon as u64);
        }
    }

    max_score = *hit_counts.0.get(&mut max_taxon).unwrap();

    while max_taxon != 0 && max_score < required_score as u64 {
        max_score = 0;
        for (taxon, count) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(max_taxon as u64, taxon as u64) {
                max_score += count;
            }
        }

        if max_score >= required_score as u64 {
            return max_taxon;
        } else {
            max_taxon = taxonomy.nodes()[max_taxon as usize].parent_id;
        }
    }

    max_taxon
}

fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }
    if &id[sz - 2..] == "/1" || &id[sz - 2..] == "/2" {
        return id[..sz - 2].to_string();
    }
    id.to_string()
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
    taxa: &mut Vec<i32>,
    hit_counts: &mut TaxonCounts,
    tx_frames: &mut Vec<i32>,
    curr_taxon_counts: &mut TaxonCounters,
) -> TaxId {
    let mut call = 0;
    taxa.clear();
    hit_counts.0.clear();
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
            if opts.use_translated_search {
                scanner.load_sequence(&tx_frames[frame_idx].to_string(), 0, 5);
            } else {
                scanner.load_sequence(if mate_num == 0 { &dna.seq } else { &dna2.seq }, 0, 5);
            }
            let mut last_minimizer = u64::MAX;
            let mut last_taxon = TAXID_MAX;

            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon = if scanner.is_ambiguous() {
                    AMBIGUOUS_SPAN_TAXON
                } else {
                    if *minimizer != last_minimizer {
                        let mut skip_lookup = false;
                        if let Some(minimum_acceptable_hash_value) =
                            idx_opts.minimum_acceptable_hash_value
                        {
                            if kv_store::murmurhash3(*minimizer) < minimum_acceptable_hash_value {
                                skip_lookup = true;
                            }
                        }

                        let taxon = if !skip_lookup {
                            hash.get(*minimizer)
                        } else {
                            Some(0)
                        };

                        last_taxon = taxon;
                        last_minimizer = *minimizer;

                        if taxon != 0 {
                            minimizer_hit_groups += 1;
                            curr_taxon_counts
                                .get_mut(&taxon)
                                .unwrap()
                                .add_kmer(scanner.next_minimizer());
                        }

                        taxon
                    } else {
                        Some(last_taxon) as u64
                    }
                };

                if taxon != 0 {
                    if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                        call = taxon;
                        break;
                    }

                    *hit_counts.0.entry(taxon).or_insert(0) += 1;
                }

                taxa.push(taxon);
            }

            if opts.use_translated_search && frame_idx != 5 {
                taxa.push(READING_FRAME_BORDER_TAXON as i32);
            }
        }

        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON as i32);
        }

        if call != 0 {
            break;
        }
    }

    let mut total_kmers = taxa.len();
    if opts.paired_end_processing {
        total_kmers -= 1;
    }
    if opts.use_translated_search {
        total_kmers -= if opts.paired_end_processing { 4 } else { 2 };
    }
    call = resolve_tree(hit_counts, taxonomy, total_kmers, opts);

    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    if call != 0 {
        stats.total_classified += 1;
        curr_taxon_counts
            .get_mut(&call)
            .unwrap()
            .increment_read_count();
    }

    if call {
        // write!(koss, "C\t").unwrap();
        koss.push_str("C\t")
    } else {
        // write!(koss, "U\t").unwrap();
        koss.push_str("U\t")
    }
    if !opts.paired_end_processing {
        koss.push_str(&format!("{}\t", dna.id));
        // write!(koss, "{}\t", dna.id).unwrap();
    } else {
        // write!(koss, "{}\t", trim_pair_info(&dna.id)).unwrap();
        koss.push_str(&format!("{}\t", trim_pair_info(&dna.id)));
    }

    let ext_call = taxonomy.nodes()[call as usize].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            Some(&taxonomy.name_data()[taxonomy.nodes()[call as usize].name_offset..])
        } else {
            None
        };

        // write!(
        //     &mut koss,
        //     "{} (taxid {})",
        //     name.unwrap_or("unclassified"),
        //     ext_call
        // )
        koss.push_str(format!("{} (taxid {})", name.unwrap_or("unclassified"), ext_call).as_str());
    } else {
        // write!(koss, "{}", ext_call).unwrap();
        koss.push_str(format!("{}", ext_call).as_str());
    }

    // write!(koss, "\t").unwrap();
    koss.push_str("\t");

    if !opts.paired_end_processing {
        koss.push_str(&format!("{}\t", dna.seq.len()));
    } else {
        koss.push_str(&format!("{}|{}\t", dna.seq.len(), dna2.seq.len()));
    }

    if opts.quick_mode {
        koss.push_str(&format!("{}:Q", ext_call));
    } else {
        if taxa.is_empty() {
            koss.push_str("0:0");
        } else {
            let mut taxa_u64: Vec<u64> = taxa.iter().map(|&x| x as u64).collect();
            add_hitlist_string(&mut koss, &mut taxa_u64, taxonomy);
        }
    }
    koss.push_str("\n");
    // write!(koss, "\n").unwrap();

    call
}

fn add_hitlist_string(oss: &mut String, taxa: &mut Vec<TaxId>, taxonomy: &Taxonomy) {
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in &taxa[1..] {
        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    write!(oss, "A:{} ", code_count).unwrap();
                } else {
                    let ext_code = taxonomy.nodes()[last_code as usize].external_id;
                    write!(oss, "{}:{} ", ext_code, code_count).unwrap();
                }
            } else {
                write!(
                    oss,
                    "{} ",
                    if last_code == MATE_PAIR_BORDER_TAXON {
                        "|:|"
                    } else {
                        "-:-"
                    }
                )
                .unwrap();
            }
            code_count = 1;
            last_code = code;
        }
    }
    if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
        if last_code == AMBIGUOUS_SPAN_TAXON {
            write!(oss, "A:{} ", code_count).unwrap();
        } else {
            let ext_code = taxonomy.nodes()[last_code as usize].external_id;
            write!(oss, "{}:{}", ext_code, code_count).unwrap();
        }
    } else {
        write!(
            oss,
            "{}",
            if last_code == MATE_PAIR_BORDER_TAXON {
                "|:|"
            } else {
                "-:-"
            }
        )
        .unwrap();
    }
}

fn initialize_outputs(opts: &Options, outputs: &mut OutputStreamData) {
    if !outputs.initialized {
        if !opts.classified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<&str> = opts.classified_output_filename.split('#').collect();
                if fields.len() < 2 {
                    panic!(
                        "Paired filename format missing # character: {}",
                        opts.classified_output_filename
                    );
                } else if fields.len() > 2 {
                    panic!(
                        "Paired filename format has >1 # character: {}",
                        opts.classified_output_filename
                    );
                }
                outputs.classified_output1 =
                    Some(File::create(format!("{}_1{}", fields[0], fields[1])).unwrap());
                outputs.classified_output2 =
                    Some(File::create(format!("{}_2{}", fields[0], fields[1])).unwrap());
            } else {
                outputs.classified_output1 =
                    Some(File::create(&opts.classified_output_filename).unwrap());
            }
            outputs.printing_sequences = true;
        }
        if !opts.unclassified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<&str> = opts.unclassified_output_filename.split('#').collect();
                if fields.len() < 2 {
                    panic!(
                        "Paired filename format missing # character: {}",
                        opts.unclassified_output_filename
                    );
                } else if fields.len() > 2 {
                    panic!(
                        "Paired filename format has >1 # character: {}",
                        opts.unclassified_output_filename
                    );
                }
                outputs.unclassified_output1 =
                    Some(File::create(format!("{}_1{}", fields[0], fields[1])).unwrap());
                outputs.unclassified_output2 =
                    Some(File::create(format!("{}_2{}", fields[0], fields[1])).unwrap());
            } else {
                outputs.unclassified_output1 =
                    Some(File::create(&opts.unclassified_output_filename).unwrap());
            }
            outputs.printing_sequences = true;
        }
        if !opts.kraken_output_filename.is_empty() {
            if opts.kraken_output_filename == Some("-") {
                outputs.kraken_output = None;
            } else {
                outputs.kraken_output = Some(File::create(&opts.kraken_output_filename).unwrap());
            }
        }
        outputs.initialized = true;
    }
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: i32) {
    match dna.format {
        SequenceFormat::Fastq => {
            if dna.seq.len() != dna.quals.len() {
                panic!(
                    "{}: Sequence length ({}) != Quality string length ({})",
                    dna.id,
                    dna.seq.len(),
                    dna.quals.len()
                );
            }
        }
        _ => return,
    }
    for i in 0..dna.seq.len() {
        if (dna.quals[i] as i32 - '!' as i32) < minimum_quality_score {
            dna.seq[i] = 'x';
        }
    }
}

fn parse_command_line(opts: &mut Options) {
    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => {
                usage(0);
            }
            "-H" => {
                i += 1;
                opts.index_filename = args[i].clone();
            }
            "-t" => {
                i += 1;
                opts.taxonomy_filename = args[i].clone();
            }
            "-T" => {
                i += 1;
                opts.confidence_threshold = args[i].parse().unwrap();
                if opts.confidence_threshold < 0.0 || opts.confidence_threshold > 1.0 {
                    panic!("confidence threshold must be in [0, 1]");
                }
            }
            "-o" => {
                i += 1;
                opts.options_filename = args[i].clone();
            }
            "-q" => {
                opts.quick_mode = true;
            }
            "-p" => {
                i += 1;
                opts.num_threads = args[i].parse().unwrap();
                if opts.num_threads < 1 {
                    panic!("number of threads can't be less than 1");
                }
            }
            "-g" => {
                i += 1;
                opts.minimum_hit_groups = args[i].parse().unwrap();
            }
            "-P" => {
                opts.paired_end_processing = true;
            }
            "-S" => {
                opts.paired_end_processing = true;
                opts.single_file_pairs = true;
            }
            "-m" => {
                opts.mpa_style_report = true;
            }
            "-K" => {
                opts.report_kmer_data = true;
            }
            "-R" => {
                i += 1;
                opts.report_filename = Some(args[i].clone());
            }
            "-z" => {
                opts.report_zero_counts = true;
            }
            "-C" => {
                i += 1;
                opts.classified_output_filename = Some(args[i].clone());
            }
            "-U" => {
                i += 1;
                opts.unclassified_output_filename = Some(args[i].clone());
            }
            "-O" => {
                i += 1;
                opts.kraken_output_filename = Some(args[i].clone());
            }
            "-n" => {
                opts.print_scientific_name = true;
            }
            "-Q" => {
                i += 1;
                opts.minimum_quality_score = args[i].parse().unwrap();
            }
            "-M" => {
                opts.use_memory_mapping = true;
            }
            _ => {}
        }
        i += 1;
    }

    if opts.index_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        eprintln!("mandatory filename missing");
        usage(0);
    }

    if opts.mpa_style_report && opts.report_filename.is_empty() {
        eprintln!("-m requires -R be used");
        usage(0);
    }
}

fn usage(exit_code: i32) {
    writeln!(
        io::stderr(),
        "Usage: classify [options] <fasta/fastq file(s)>\n
    Options: (*mandatory)\n
    * -H filename      Kraken 2 index filename\n
    * -t filename      Kraken 2 taxonomy filename\n
    * -o filename      Kraken 2 options filename\n
      -q               Quick mode\n
      -M               Use memory mapping to access hash & taxonomy\n
      -T NUM           Confidence score threshold (def. 0)\n
      -p NUM           Number of threads (def. 1)\n
      -Q NUM           Minimum quality score (FASTQ only, def. 0)\n
      -P               Process pairs of reads\n
      -S               Process pairs with mates in same file\n
      -R filename      Print report to filename\n
      -m               In comb. w/ -R, use mpa-style report\n
      -z               In comb. w/ -R, report taxa w/ 0 count\n
      -n               Print scientific name instead of taxid in Kraken output\n
      -g NUM           Minimum number of hit groups needed for call\n
      -C filename      Filename/format to have classified sequences\n
      -U filename      Filename/format to have unclassified sequences\n
      -O filename      Output file for normal Kraken output\n
      -K               In comb. w/ -R, provide minimizer information in report"
    )
    .unwrap();
    process::exit(exit_code);
}

fn update_total_counters(
    total_taxon_counters: Arc<Mutex<TaxonCounters>>,
    thread_taxon_counters: TaxonCounters,
) {
    let mut total_counters = total_taxon_counters.lock().unwrap();
    for (key, value) in thread_taxon_counters {
        *total_counters.entry(key).or_insert(0) += value;
    }
}

fn output_results(
    outputs: Arc<Mutex<BinaryHeap<OutputData>>>,
    output_file_path: &str,
) -> io::Result<()> {
    let outputs = outputs.lock().unwrap();
    let mut output_file = File::create(output_file_path)?;

    for output_data in outputs.iter() {
        writeln!(output_file, "{}", output_data.kraken_str)?;
        // Write other output data fields as needed
    }

    Ok(())
}

fn read_input<P: AsRef<Path>>(
    filename1: Option<P>,
    filename2: Option<P>,
    // Assume options determine how inputs are processed, e.g., paired-end
) -> io::Result<(
    SequenceReader<BufReader<File>>,
    Option<SequenceReader<BufReader<File>>>,
)> {
    let fptr1 = match filename1 {
        Some(f) => BufReader::new(File::open(f)?),
        None => BufReader::new(io::stdin()),
    };

    let reader1 = SequenceReader::new(fptr1);

    let reader2 = if let Some(f) = filename2 {
        let fptr2 = BufReader::new(File::open(f)?);
        Some(SequenceReader::new(fptr2))
    } else {
        None
    };

    Ok((reader1, reader2))
}

fn process_inputs<P: AsRef<Path>>(
    filename1: Option<P>,
    filename2: Option<P>,
    options: Options,
) -> io::Result<()> {
    // Setup input readers
    let (mut reader1, reader2) = read_input(filename1, filename2)?;

    // Setup shared state for statistics (assuming concurrent processing in a real scenario)
    let classification_stats = Arc::new(Mutex::new(ClassificationStats {
        total_sequences: 0,
        total_bases: todo!(),
        total_classified: todo!(),
    }));

    // If paired-end processing is enabled and a second reader is available, use it
    let maybe_reader2 = if options.paired_end_processing {
        reader2.as_mut()
    } else {
        None
    };

    // Process fragments
    process_fragments(&mut reader1, maybe_reader2, classification_stats.clone())?;

    // In a real application, here we would likely perform additional steps, such as:
    // - Updating global counters or results based on `classification_stats`
    // - Writing output to files or stdout
    // - Handling any cleanup or finalization tasks

    Ok(())
}

fn open_file(filename: Option<&str>) -> io::Result<Box<dyn BufRead>> {
    match filename {
        Some(path) => Ok(Box::new(BufReader::new(File::open(path)?))),
        None => Ok(Box::new(BufReader::new(io::stdin()))),
    }
}
// Implement `Ord` and `PartialOrd` to establish how `OutputData` instances are compared.
impl PartialEq for OutputData {
    fn eq(&self, other: &Self) -> bool {
        // Implement the equality comparison logic here
        // Return true if the two instances are equal, false otherwise
        todo!()
    }
}

impl Eq for OutputData {}

impl Ord for OutputData {
    fn cmp(&self, other: &Self) -> Ordering {
        todo!()
    }

    fn max(self, other: Self) -> Self
    where
        Self: Sized,
    {
        std::cmp::max_by(self, other, Ord::cmp)
    }

    fn min(self, other: Self) -> Self
    where
        Self: Sized,
    {
        std::cmp::min_by(self, other, Ord::cmp)
    }

    fn clamp(self, min: Self, max: Self) -> Self
    where
        Self: Sized,
        Self: PartialOrd,
    {
        assert!(min <= max);
        if self < std::cmp::min {
            std::cmp::min
        } else if self > std::cmp::max {
            std::cmp::max
        } else {
            self
        }
    }
    // Rest of the code...
}

impl PartialOrd for OutputData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn initialize_output_queue() -> BinaryHeap<OutputData> {
    BinaryHeap::new()
}

/*
#[derive(Eq, PartialEq)]
Remove the duplicate definition of OutputData
struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}
*/

// impl Ord for OutputData {
//     fn cmp(&self, other: &Self) -> Ordering {
//         // This comparison is reversed to turn the BinaryHeap (max heap by default) into a min heap.
//         other.block_id.cmp(&self.block_id)
//     }
// }

// impl PartialOrd for OutputData {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         Some(self.cmp(other))
//     }
// }

// impl OutputData {
//     // Constructor for creating a new instance of OutputData
//     pub fn new(
//         block_id: u64,
//         kraken_str: String,
//         classified_out1_str: String,
//         classified_out2_str: String,
//         unclassified_out1_str: String,
//         unclassified_out2_str: String,
//     ) -> Self {
//         OutputData {
//             block_id,
//             kraken_str,
//             classified_out1_str,
//             classified_out2_str,
//             unclassified_out1_str,
//             unclassified_out2_str,
//         }
//     }
// }
fn process_files<P: AsRef<Path>>(
    filename1: Option<P>,
    filename2: Option<P>,
    options: Options,
) -> io::Result<()> {
    // Open input files or stdin for reading
    let reader1 = open_file(filename1.as_ref().map(|p| p.as_ref().to_str().unwrap()))?;
    let reader2 = if options.paired_end_processing {
        filename2.map(|p| open_file(Some(p.as_ref().to_str().unwrap())).unwrap())
    } else {
        None
    };

    // Initialize the output queue and other necessary data structures
    let output_queue = Arc::new(Mutex::new(initialize_output_queue()));
    let classification_stats = Arc::new(Mutex::new(ClassificationStats {
        total_sequences: 0,
        total_bases: todo!(),
        total_classified: todo!(),
    }));

    let total_taxon_counters = Arc::new(Mutex::new(TaxonCounters::new()));

    // Process the input fragments
    let maybe_reader2 = reader2.as_ref().map(|r| r.as_ref());
    process_fragments(&mut *reader1, maybe_reader2, classification_stats.clone())?;

    // Update total counters based on the results from `process_fragments`
    // Assume this function takes care of updating `total_taxon_counters`
    // based on local counters from processing
    // update_total_counters(total_taxon_counters, ...);

    // Output results
    // Assume output paths or writers are defined elsewhere
    output_results(output_queue, "classified_out1_path")?;

    Ok(())
}

// fn main() -> io::Result<()> {
//     let options = Options {
//         paired_end_processing: true,
//         // Initialize other options as needed
//     };

//     // Call process_files with example file paths or None for stdin
//     Ok(process_files(
//         Some(Path::new("path/to/file1")),
//         Some(Path::new("path/to/file2")),
//         options,
//     ))
// }

// Assuming a function for sequence classification exists
// fn classify_sequence(seq: &Sequence) -> bool { /* Dummy implementation */ true }

// Simplified SequenceReader for demonstration
struct SequenceReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> SequenceReader<R> {
    fn new(reader: R) -> Self {
        Self { reader }
    }

    fn next_sequence(&mut self) -> Result<Option<Sequence>, io::Error> {
        // Simplified: Pretend to read and return the next sequence or None if at EOF
        Ok(None) // Placeholder implementation
    }
}

fn process_fragments<R: BufRead>(
    reader1: &mut SequenceReader<R>,
    reader2: Option<&mut SequenceReader<R>>,
    classification_stats: Arc<Mutex<ClassificationStats>>,
) -> Result<(), io::Error> {
    loop {
        let seq1 = match reader1.next_sequence()? {
            Some(seq) => seq,
            None => break, // EOF reached for reader1
        };

        let seq2 = if let Some(ref mut r2) = reader2 {
            match r2.next_sequence()? {
                Some(seq) => Some(seq),
                None => None, // EOF reached or not paired-end processing
            }
        } else {
            None
        };

        // Example classification logic (simplified)
        let classified = classify_sequence(
            dna,
            dna2,
            koss,
            hash,
            taxonomy,
            idx_opts,
            opts,
            stats,
            scanner,
            taxa,
            hit_counts,
            tx_frames,
            curr_taxon_counts,
        );
        {
            let mut stats = classification_stats.lock().unwrap();
            stats.total_sequences += 1;
            if classified {
                stats.total_classified += 1;
            }
        }

        // Further processing based on classification...
        // Update output data, etc.
    }

    Ok(())
}

// Placeholder function signatures for the defined operations:
// fn open_file(filename: Option<&str>) -> io::Result<Box<dyn BufRead>>;
// fn initialize_output_queue() -> BinaryHeap<OutputData>;
// fn read_input<P: AsRef<Path>>(...);
// fn process_fragments<R: BufRead>(...);
// fn update_total_counters(...);
// fn output_results(...);

// Assuming these structs are defined according to previous instructions

fn process_files<P: AsRef<Path>>(
    filename1: Option<P>,
    filename2: Option<P>,
    options: Options,
) -> io::Result<()> {
    // Open input files or stdin for reading
    let reader1 = open_file(filename1.as_ref().map(|p| p.as_ref().to_str().unwrap()))?;
    let reader2 = if options.paired_end_processing {
        filename2.map(|p| open_file(Some(p.as_ref().to_str().unwrap())).unwrap())
    } else {
        None
    };

    // Initialize the output queue and other necessary data structures
    let output_queue = Arc::new(Mutex::new(initialize_output_queue()));
    let classification_stats = Arc::new(Mutex::new(ClassificationStats { total_sequences: 0 }));
    let total_taxon_counters = Arc::new(Mutex::new(TaxonCounters::new()));

    // Process the input fragments
    let maybe_reader2 = reader2.as_ref().map(|r| r.as_ref());
    process_fragments(&mut *reader1, maybe_reader2, classification_stats.clone())?;

    // Update total counters based on the results from `process_fragments`
    // Assume this function takes care of updating `total_taxon_counters`
    // based on local counters from processing
    // update_total_counters(total_taxon_counters, ...);

    // Output results
    // Assume output paths or writers are defined elsewhere
    output_results(
        output_queue,
        "classified_out1_path",
        "classified_out2_path",
        "unclassified_out1_path",
        "unclassified_out2_path",
    )?;

    Ok(())
}
