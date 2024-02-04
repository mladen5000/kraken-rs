use crate::aa_translate::translate_to_all_frames;
use crate::kv_store::KeyValueStore;
// use crate::build_db::TaxId;
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::IndexOptions;
use crate::kraken2_data::TaxId;
use crate::kraken2_data::TAXID_MAX;
use crate::kv_store;
use crate::kv_store::*;
use crate::reports;
use crate::reports::{TaxonCounters, TaxonCounts};
use num_cpus;
use std::process;

use crate::seqreader::{Sequence, SequenceFormat};
use crate::taxonomy::*;
use std::collections::{BinaryHeap, HashMap};
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::io::{BufReader, Read};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

pub const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
pub const MATE_PAIR_BORDER_TAXON: TaxId = TAXID_MAX;
pub const READING_FRAME_BORDER_TAXON: TaxId = TAXID_MAX - 1;
pub const AMBIGUOUS_SPAN_TAXON: TaxId = TAXID_MAX - 2;

pub struct Options {
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
pub struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}
pub struct OutputStreamData {
    classified_output1: Option<File>,
    classified_output2: Option<File>,
    unclassified_output1: Option<File>,
    unclassified_output2: Option<File>,
    kraken_output: Option<File>,
    initialized: bool,
    printing_sequences: bool,
}
pub struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}

fn main() {
    let mut opts = Options::default();
    let mut taxon_counters = TaxonCounters::new(); // stats per taxon

    parse_command_line(&mut opts);

    println!("Loading database information...");

    let mut idx_opts = IndexOptions::default();
    let opts_filesize = std::fs::metadata(&opts.options_filename).unwrap().len();
    let mut idx_opt_fs = std::fs::File::open(&opts.options_filename).unwrap();
    unsafe {
        std::ptr::copy_nonoverlapping(
            &mut idx_opts as *mut _ as *mut _,
            idx_opt_fs.as_mut_ptr() as *mut _,
            opts_filesize as usize,
        );
    }
    opts.use_translated_search = !idx_opts.dna_db;

    let mut taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping);
    let mut hash_ptr = CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping);

    println!(" done.");

    let mut stats = ClassificationStats::default();

    let mut outputs = OutputStreamData::default();

    let start_time = Instant::now();
    if std::env::args().len() == 1 {
        if opts.paired_end_processing && !opts.single_file_pairs {
            panic!("paired end processing used with no files specified");
        }
        process_files(
            None,
            None,
            &mut hash_ptr,
            &mut taxonomy,
            &mut idx_opts,
            &mut opts,
            &mut stats,
            &mut outputs,
            &mut taxon_counters,
        );
    } else {
        for i in 1..std::env::args().len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == std::env::args().len() {
                    panic!("paired end processing used with unpaired file");
                }
                process_files(
                    Some(std::env::args().nth(i).unwrap()),
                    Some(std::env::args().nth(i + 1).unwrap()),
                    &mut hash_ptr,
                    &mut taxonomy,
                    &mut idx_opts,
                    &mut opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                );
            } else {
                process_files(
                    Some(std::env::args().nth(i).unwrap()),
                    None,
                    &mut hash_ptr,
                    &mut taxonomy,
                    &mut idx_opts,
                    &mut opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                );
            }
        }
    }

    let end_time = Instant::now();

    report_stats(start_time, end_time, &mut stats);

    if !opts.report_filename.is_empty() {
        if opts.mpa_style_report {
            reports::report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &mut taxonomy,
                &mut taxon_counters,
            );
        } else {
            let total_unclassified = stats.total_sequences - stats.total_classified;
            reports::report_kraken_style(
                &opts.report_filename,
                opts.report_zero_counts,
                opts.report_kmer_data,
                &mut taxonomy,
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

fn process_files(
    filename1: &str,
    filename2: &str,
    hash: &mut KeyValueStore,
    tax: &mut Taxonomy,
    idx_opts: &mut IndexOptions,
    opts: &mut Options,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    total_taxon_counters: &mut TaxonCounters,
) {
    let fptr1 = match filename1 {
        Some(name) => Box::new(BufReader::new(File::open(name).unwrap())) as Box<dyn Read>,
        None => Box::new(BufReader::new(std::io::stdin())) as Box<dyn Read>,
    };

    let fptr2 = match filename2 {
        Some(name) => Some(Box::new(BufReader::new(File::open(name).unwrap())) as Box<dyn Read>),
        None => None,
    };

    let output_queue: Arc<Mutex<BinaryHeap<OutputData>>> = Arc::new(Mutex::new(BinaryHeap::new()));

    let threads: Vec<_> = (0..num_cpus::get())
        .map({
            let queue_clone = Arc::clone(&output_queue);
            let fptr1_clone = fptr1.clone();
            let fptr2_clone = fptr2.clone();

            thread::spawn(move || {
                // Thread-local variables go here

                loop {
                    // Read and process a block of data from fptr1 and possibly fptr2
                    let mut thread_stats = ClassificationStats {
                        total_sequences: 0,
                        total_bases: 0,
                        total_classified: 0,
                    };
                    let mut kraken_oss = String::new();
                    let mut c1_oss = String::new();
                    let mut c2_oss = String::new();
                    let mut u1_oss = String::new();
                    let mut u2_oss = String::new();
                    let mut thread_taxon_counters = HashMap::new();

                    let mut ok_read = false;

                    // Input processing block
                    if !opts.paired_end_processing {
                        // Unpaired data? Just read in a sized block
                        ok_read = reader1.load_block(fptr1, 3 * 1024 * 1024);
                    } else if !opts.single_file_pairs {
                        // Paired data in 2 files? Read a line-counted batch from each file.
                        ok_read = reader1.load_batch(fptr1, NUM_FRAGMENTS_PER_THREAD);
                        if ok_read && opts.paired_end_processing {
                            ok_read = reader2.load_batch(fptr2, NUM_FRAGMENTS_PER_THREAD);
                        }
                    } else {
                        let mut frags = NUM_FRAGMENTS_PER_THREAD * 2;
                        // Ensure frag count is even - just in case above line is changed
                        if frags % 2 == 1 {
                            frags += 1;
                        }
                        ok_read = reader1.load_batch(fptr1, frags);
                    }
                    block_id = next_input_block_id.fetch_add(1, Ordering::SeqCst);

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

                    while let Some(mut seq1) = reader1.next_sequence() {
                        let mut seq2 = if opts.paired_end_processing {
                            if opts.single_file_pairs {
                                reader1.next_sequence()
                            } else {
                                reader2.next_sequence()
                            }
                        } else {
                            None
                        };

                        if seq2.is_none() {
                            break;
                        }

                        thread_stats.total_sequences += 1;

                        if opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut seq1, opts.minimum_quality_score);
                            if let Some(seq2) = &mut seq2 {
                                mask_low_quality_bases(seq2, opts.minimum_quality_score);
                            }
                        }

                        let call = classify_sequence(
                            &mut seq1,
                            &mut seq2,
                            &mut kraken_oss,
                            &mut hash,
                            &mut tax,
                            &idx_opts,
                            &opts,
                            &mut thread_stats,
                            &mut scanner,
                            &mut taxa,
                            &mut hit_counts,
                            &mut translated_frames,
                            &mut thread_taxon_counters,
                        );

                        if let Some(call) = call {
                            let buffer = format!(" kraken:taxid|{}", tax.nodes()[call].external_id);
                            seq1.header.push_str(&buffer);
                            write!(&mut c1_oss, "{}", format!("{}", seq1)).unwrap();
                            if let Some(mut seq2) = seq2 {
                                seq2.header.push_str(&buffer);
                                write!(&mut c2_oss, "{}", format!("{}", seq2)).unwrap();
                            }
                        } else {
                            write!(&mut u1_oss, "{}", format!("{}", seq1)).unwrap();
                            if let Some(seq2) = seq2 {
                                write!(&mut u2_oss, "{}", format!("{}", seq2)).unwrap();
                            }
                        }
                        thread_stats.total_bases += seq1.seq.len();
                        if let Some(seq2) = seq2 {
                            thread_stats.total_bases += seq2.seq.len();
                        }
                    }

                    stats
                        .total_sequences
                        .fetch_add(thread_stats.total_sequences, Ordering::SeqCst);
                    stats
                        .total_bases
                        .fetch_add(thread_stats.total_bases, Ordering::SeqCst);
                    stats
                        .total_classified
                        .fetch_add(thread_stats.total_classified, Ordering::SeqCst);

                    if atty::is(atty::Stream::Stderr) {
                        eprint!(
                            "\rProcessed {} sequences ({} bp) ...",
                            stats.total_sequences.load(Ordering::SeqCst),
                            stats.total_bases.load(Ordering::SeqCst)
                        );
                    }

                    if !outputs.initialized {
                        initialize_outputs(&opts, &mut outputs);
                    }

                    let mut out_data = OutputData {
                        block_id: block_id,
                        kraken_str: kraken_oss.clone(),
                        classified_out1_str: c1_oss.clone(),
                        classified_out2_str: c2_oss.clone(),
                        unclassified_out1_str: u1_oss.clone(),
                        unclassified_out2_str: u2_oss.clone(),
                    };

                    {
                        let mut output_queue_guard = output_queue.lock().unwrap();
                        output_queue_guard.push(out_data);
                    }

                    {
                        let mut total_taxon_counters_guard = total_taxon_counters.lock().unwrap();
                        for (key, value) in thread_taxon_counters.drain() {
                            *total_taxon_counters_guard.entry(key).or_insert(0) += value;
                        }
                    }

                    let mut output_loop = true;

                    while output_loop {
                        let mut next_output_block_id_guard = next_output_block_id.lock().unwrap();
                        let mut output_queue_guard = output_queue.lock().unwrap();
                        output_loop = !output_queue_guard.is_empty();
                        if output_loop {
                            out_data = output_queue_guard.peek().unwrap().clone();
                            if out_data.block_id == *next_output_block_id_guard {
                                output_queue_guard.pop();
                                *next_output_block_id_guard += 1;
                            } else {
                                output_loop = false;
                            }
                        }
                        drop(next_output_block_id_guard);
                        drop(output_queue_guard);

                        if !output_loop {
                            break;
                        }

                        if let Some(kraken_output) = &mut outputs.kraken_output {
                            write!(kraken_output, "{}", out_data.kraken_str).unwrap();
                        }
                        if let Some(classified_output1) = &mut outputs.classified_output1 {
                            write!(classified_output1, "{}", out_data.classified_out1_str).unwrap();
                        }
                        if let Some(classified_output2) = &mut outputs.classified_output2 {
                            write!(classified_output2, "{}", out_data.classified_out2_str).unwrap();
                        }
                        if let Some(unclassified_output1) = &mut outputs.unclassified_output1 {
                            write!(unclassified_output1, "{}", out_data.unclassified_out1_str)
                                .unwrap();
                        }
                        if let Some(unclassified_output2) = &mut outputs.unclassified_output2 {
                            write!(unclassified_output2, "{}", out_data.unclassified_out2_str)
                                .unwrap();
                        }
                    }
                }
            })
        })
        .collect();

    for t in threads {
        t.join().unwrap();
    }

    // Flush the output streams
}

fn resolve_tree(
    hit_counts: &mut TaxonCounts,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (opts.confidence_threshold * total_minimizers as f32).ceil() as u32;

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

    max_score = hit_counts[&max_taxon];

    while max_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(max_taxon, taxon) {
                max_score += count;
            }
        }

        if max_score >= required_score {
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
    taxa: &mut Vec<_>,
    hit_counts: &mut TaxonCounts,
    tx_frames: &mut Vec<_>,
    curr_taxon_counts: &mut TaxonCounters,
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

            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon = if scanner.is_ambiguous() {
                    AMBIGUOUS_SPAN_TAXON
                } else {
                    if minimizer != last_minimizer {
                        let mut skip_lookup = false;
                        if let Some(minimum_acceptable_hash_value) =
                            idx_opts.minimum_acceptable_hash_value
                        {
                            if kv_store::murmurhash3(minimizer) < minimum_acceptable_hash_value {
                                skip_lookup = true;
                            }
                        }

                        let taxon = if !skip_lookup { hash.get(minimizer) } else { 0 };

                        last_taxon = taxon;
                        last_minimizer = minimizer;

                        if taxon != 0 {
                            minimizer_hit_groups += 1;
                            curr_taxon_counts
                                .get_mut(&taxon)
                                .unwrap()
                                .add_kmer(scanner.last_minimizer());
                        }

                        taxon
                    } else {
                        last_taxon
                    }
                };

                if taxon != 0 {
                    if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                        call = taxon;
                        break;
                    }

                    *hit_counts.entry(taxon).or_insert(0) += 1;
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
        write!(koss, "C\t").unwrap();
    } else {
        write!(koss, "U\t").unwrap();
    }
    if !opts.paired_end_processing {
        write!(koss, "{}\t", dna.id).unwrap();
    } else {
        write!(koss, "{}\t", trim_pair_info(&dna.id)).unwrap();
    }

    let ext_call = taxonomy.nodes()[call as usize].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            Some(&taxonomy.name_data()[taxonomy.nodes()[call as usize].name_offset..])
        } else {
            None
        };

        write!(
            koss,
            "{} (taxid {})",
            name.unwrap_or("unclassified"),
            ext_call
        )
        .unwrap();
    } else {
        write!(koss, "{}", ext_call).unwrap();
    }

    write!(koss, "\t").unwrap();

    if !opts.paired_end_processing {
        write!(koss, "{}\t", dna.seq.len()).unwrap();
    } else {
        write!(koss, "{}|{}\t", dna.seq.len(), dna2.seq.len()).unwrap();
    }

    if opts.quick_mode {
        write!(koss, "{}:Q", ext_call).unwrap();
    } else {
        if taxa.is_empty() {
            write!(koss, "0:0").unwrap();
        } else {
            add_hitlist_string(koss, taxa, taxonomy);
        }
    }

    write!(koss, "\n").unwrap();

    call
}

fn add_hitlist_string(oss: &mut String, taxa: &Vec<TaxId>, taxonomy: &Taxonomy) {
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
            if opts.kraken_output_filename == "-" {
                outputs.kraken_output = None;
            } else {
                outputs.kraken_output = Some(File::create(&opts.kraken_output_filename).unwrap());
            }
        }
        outputs.initialized = true;
    }
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: i32) {
    if dna.format != Format::Fastq {
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
                opts.report_filename = args[i].clone();
            }
            "-z" => {
                opts.report_zero_counts = true;
            }
            "-C" => {
                i += 1;
                opts.classified_output_filename = args[i].clone();
            }
            "-U" => {
                i += 1;
                opts.unclassified_output_filename = args[i].clone();
            }
            "-O" => {
                i += 1;
                opts.kraken_output_filename = args[i].clone();
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
