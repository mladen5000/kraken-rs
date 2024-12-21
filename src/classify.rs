// lib.rs or main.rs

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Write};
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

// Placeholder types, assuming these are defined elsewhere
type TaxonCountsT = HashMap<taxid_t, u32>;
type TaxonCountersT = HashMap<taxid_t, readcounts::TaxonCount>; // or appropriate struct

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

pub fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: &KeyValueStore,
    tax: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &Arc<Mutex<ClassificationStats>>,
    outputs: &Arc<Mutex<OutputStreamData>>,
    total_taxon_counters: &Arc<Mutex<TaxonCountersT>>,
) {
    let fptr1: Box<dyn BufRead> = if let Some(fname) = filename1 {
        Box::new(BufReader::new(
            File::open(fname).expect("Unable to open file1"),
        ))
    } else {
        Box::new(BufReader::new(stdin()))
    };

    let fptr2: Option<Box<dyn BufRead>> = if opts.paired_end_processing && !opts.single_file_pairs {
        Some(Box::new(BufReader::new(
            File::open(filename2.expect("No paired file provided")).expect("Unable to open file2"),
        )))
    } else {
        None
    };

    // We will create channels for communication:
    // - A channel from main reading thread to workers: (block_id, batch_of_sequences)
    // - A channel from workers back to main output thread: OutputData
    let (work_tx, work_rx) = unbounded::<(u64, Vec<(Sequence, Option<Sequence>)>)>();
    let (result_tx, result_rx) = unbounded::<OutputData>();

    // Reading thread
    // This thread is responsible for reading sequences in blocks and sending them to workers.
    let opts_clone = opts.clone();
    let paired = opts.paired_end_processing;
    let single_file_pairs = opts.single_file_pairs;
    let mut block_id_counter = 0u64;

    let handle_reader = {
        let opts_clone = opts_clone;
        let mut reader1 = BatchequenceReader::new();
        let mut reader2 = if paired && !single_file_pairs {
            Some(BatchequenceReader::new())
        } else {
            None
        };

        let mut fptr1 = fptr1;
        let mut fptr2 = fptr2;
        thread::spawn(move || {
            loop {
                let ok_read = if !paired {
                    // Unpaired data: load a chunk of ~3MB
                    reader1.load_block(&mut *fptr1, 3 * 1024 * 1024)
                } else if !single_file_pairs {
                    // Paired data in 2 files: load fixed # of fragments from each
                    let r1 = reader1.load_batch(&mut *fptr1, NUM_FRAGMENTS_PER_THREAD);
                    let r2 = if r1 {
                        reader2
                            .as_mut()
                            .unwrap()
                            .load_batch(&mut *fptr2.as_mut().unwrap(), NUM_FRAGMENTS_PER_THREAD)
                    } else {
                        false
                    };
                    r1 && r2
                } else {
                    // Paired data in single file: load double number of sequences
                    let frags = NUM_FRAGMENTS_PER_THREAD * 2;
                    reader1.load_batch(&mut *fptr1, frags)
                };

                if !ok_read {
                    // no more data
                    break;
                }

                let mut sequences = Vec::new();
                loop {
                    let mut seq1 = Sequence::default();
                    let valid_fragment = reader1.next_sequence(&mut seq1);
                    if paired && valid_fragment {
                        let mut seq2 = Sequence::default();
                        let valid_fragment = if single_file_pairs {
                            reader1.next_sequence(&mut seq2)
                        } else {
                            reader2.as_mut().unwrap().next_sequence(&mut seq2)
                        };
                        if !valid_fragment {
                            break;
                        }
                        sequences.push((seq1, Some(seq2)));
                    } else if valid_fragment {
                        sequences.push((seq1, None));
                    } else {
                        break;
                    }
                }
                let current_block_id = block_id_counter;
                block_id_counter += 1;
                work_tx.send((current_block_id, sequences)).unwrap();
            }
            // Signal no more work
            drop(work_tx);
        })
    };

    // Worker threads: use a thread pool. For simplicity, spawn a fixed number of threads equal to opts.num_threads.
    // Each worker receives blocks of sequences, classifies them, and sends results back.
    let taxonomy = Arc::new(tax.clone());
    let hash_arc = Arc::new(hash.clone());
    let idx_opts_arc = Arc::new(idx_opts.clone());
    let opts_arc = Arc::new(opts.clone());

    let mut handles = Vec::new();
    for _ in 0..opts.num_threads {
        let work_rx = work_rx.clone();
        let result_tx = result_tx.clone();
        let taxonomy = Arc::clone(&taxonomy);
        let hash_arc = Arc::clone(&hash_arc);
        let idx_opts_arc = Arc::clone(&idx_opts_arc);
        let opts_arc = Arc::clone(&opts_arc);
        let stats_arc = Arc::clone(&stats);
        let total_taxon_counters = Arc::clone(&total_taxon_counters);
        let outputs_arc = Arc::clone(&outputs);

        handles.push(thread::spawn(move || {
            let mut scanner = MinimizerScanner::new(
                idx_opts_arc.k,
                idx_opts_arc.l,
                idx_opts_arc.spaced_seed_mask,
                idx_opts_arc.dna_db,
                idx_opts_arc.toggle_mask,
                idx_opts_arc.revcom_version,
            );
            let mut taxa = Vec::new();
            let mut hit_counts: TaxonCountsT = TaxonCountsT::default();
            let mut translated_frames = vec![String::new(); 6];
            let mut local_taxon_counters: TaxonCountersT = TaxonCountersT::default();

            loop {
                let (block_id, sequences) = match work_rx.recv() {
                    Ok(data) => data,
                    Err(_) => break, // no more data
                };

                let mut thread_stats = ClassificationStats {
                    total_sequences: 0,
                    total_bases: 0,
                    total_classified: 0,
                };

                let mut kraken_str = String::new();
                let mut classified_out1_str = String::new();
                let mut classified_out2_str = String::new();
                let mut unclassified_out1_str = String::new();
                let mut unclassified_out2_str = String::new();

                // Lazy initialization of outputs, done once
                {
                    let mut outputs = outputs_arc.lock().unwrap();
                    if !outputs.initialized {
                        // Just pick format from the input sequences (assuming we know it)
                        let format = sequences
                            .first()
                            .map(|(s, _)| s.format)
                            .unwrap_or(SequenceFormat::Fasta);
                        InitializeOutputs(&opts_arc, &mut outputs, format);
                    }
                }

                for (mut seq1, seq2_opt) in sequences {
                    thread_stats.total_sequences += 1;
                    if opts_arc.minimum_quality_score > 0 {
                        MaskLowQualityBases(&mut seq1, opts_arc.minimum_quality_score);
                        if let Some(ref mut seq2) = seq2_opt.as_ref().cloned() {
                            let mut seq2_cloned = seq2.clone();
                            MaskLowQualityBases(&mut seq2_cloned, opts_arc.minimum_quality_score);
                        }
                    }

                    let call = ClassifySequence(
                        &mut seq1,
                        seq2_opt.as_ref().map(|s| s.clone()).unwrap_or_default(),
                        &mut kraken_str,
                        &hash_arc,
                        &taxonomy,
                        &idx_opts_arc,
                        &opts_arc,
                        &mut thread_stats,
                        &mut scanner,
                        &mut taxa,
                        &mut hit_counts,
                        &mut translated_frames,
                        &mut local_taxon_counters,
                    );

                    // Append classified/unclassified sequences
                    if call != 0 {
                        // Append to classified outputs
                        let ext_call = taxonomy.nodes()[call].external_id;
                        seq1.header.push_str(&format!(" kraken:taxid|{}", ext_call));
                        classified_out1_str.push_str(&seq1.to_string());
                        if let Some(mut seq2) = seq2_opt {
                            seq2.header.push_str(&format!(" kraken:taxid|{}", ext_call));
                            classified_out2_str.push_str(&seq2.to_string());
                        }
                    } else {
                        // Append to unclassified outputs
                        unclassified_out1_str.push_str(&seq1.to_string());
                        if let Some(seq2) = seq2_opt {
                            unclassified_out2_str.push_str(&seq2.to_string());
                        }
                    }

                    thread_stats.total_bases += seq1.seq.len() as u64;
                    if let Some(seq2) = seq2_opt {
                        thread_stats.total_bases += seq2.seq.len() as u64;
                    }

                    taxa.clear();
                    hit_counts.clear();
                }

                // Update global stats
                {
                    let mut s = stats_arc.lock().unwrap();
                    s.total_sequences += thread_stats.total_sequences;
                    s.total_bases += thread_stats.total_bases;
                    s.total_classified += thread_stats.total_classified;
                }

                // Update taxon counters
                {
                    let mut global_taxon_counters = total_taxon_counters.lock().unwrap();
                    for (k, v) in local_taxon_counters.drain() {
                        let entry = global_taxon_counters
                            .entry(k)
                            .or_insert_with(|| TaxonCounter::default());
                        entry.merge(v);
                    }
                }

                // Send output data back
                let out_data = OutputData {
                    block_id,
                    kraken_str,
                    classified_out1_str,
                    classified_out2_str,
                    unclassified_out1_str,
                    unclassified_out2_str,
                };
                result_tx.send(out_data).unwrap();
            }
        }));
    }

    drop(work_rx); // Close sending side if reading ended
    drop(result_tx); // Will close after workers done

    // Output thread: we now need to print in order
    // Since we have block_ids, we can use a min-heap or just track next_output_block_id
    let mut next_output_block_id = 0u64;
    let mut pq = BinaryHeap::new(); // We'll store Reverse(block_id) to get smallest first

    let mut outputs_locked = outputs.lock().unwrap();
    let stdout_handle = stdout();
    let kraken_output = outputs_locked.kraken_output.as_mut();
    let c_out1 = outputs_locked.classified_output1.as_mut();
    let c_out2 = outputs_locked.classified_output2.as_mut();
    let u_out1 = outputs_locked.unclassified_output1.as_mut();
    let u_out2 = outputs_locked.unclassified_output2.as_mut();
    drop(outputs_locked); // release lock

    for result in result_rx {
        // Insert into priority queue
        pq.push((Reverse(result.block_id), result));

        while let Some((Reverse(bid), out_data)) = pq.peek().cloned() {
            if bid == next_output_block_id {
                // This is the next block to print
                pq.pop();
                if let Some(kout) = kraken_output {
                    write!(kout, "{}", out_data.kraken_str).ok();
                }
                if let Some(c1) = c_out1 {
                    write!(c1, "{}", out_data.classified_out1_str).ok();
                }
                if let Some(c2) = c_out2 {
                    write!(c2, "{}", out_data.classified_out2_str).ok();
                }
                if let Some(u1) = u_out1 {
                    write!(u1, "{}", out_data.unclassified_out1_str).ok();
                }
                if let Some(u2) = u_out2 {
                    write!(u2, "{}", out_data.unclassified_out2_str).ok();
                }
                next_output_block_id += 1;
            } else {
                break;
            }
        }
    }

    // Wait for reader and workers to finish
    handle_reader.join().expect("Reader thread panicked");
    for h in handles {
        h.join().expect("Worker thread panicked");
    }

    // Flush outputs
    if let Some(kout) = kraken_output {
        kout.flush().ok();
    }
    if let Some(c1) = c_out1 {
        c1.flush().ok();
    }
    if let Some(c2) = c_out2 {
        c2.flush().ok();
    }
    if let Some(u1) = u_out1 {
        u1.flush().ok();
    }
    if let Some(u2) = u_out2 {
        u2.flush().ok();
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
