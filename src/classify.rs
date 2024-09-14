use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::TaxId;
use crate::kraken2_data::*;
use crate::kv_store::murmur_hash3;
use crate::kv_store::*;
use crate::mmscanner::*;
use crate::readcounts::*;
use crate::reports::*;
use crate::seqreader::{Sequence, SequenceFormat};
use crate::taxonomy::*;
use crate::utilities::*;

use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;

const NUM_FRAGMENTS_PER_THREAD: usize = 10_000;
const MATE_PAIR_BORDER_TAXON: TaxId = TaxId::MAX;
const READING_FRAME_BORDER_TAXON: TaxId = TaxId::MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxId = TaxId::MAX - 2;

// Corrected Type Alias
type TaxonCounters = HashMap<TaxId, TaxonCount>;

#[derive(Default, Clone)]
struct Options {
    index_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    report_filename: Option<String>,
    classified_output_filename: Option<String>,
    unclassified_output_filename: Option<String>,
    kraken_output_filename: Option<String>,
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
    minimum_quality_score: usize,
    minimum_hit_groups: usize,
    use_memory_mapping: bool,
    match_input_order: bool,
}

#[derive(Default, Clone)]
struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

#[derive(Default)]
struct OutputStreamData {
    initialized: bool,
    kraken_output: Option<BufWriter<Box<dyn Write + Send>>>,
    classified_output1: Option<BufWriter<Box<dyn Write + Send>>>,
    classified_output2: Option<BufWriter<Box<dyn Write + Send>>>,
    unclassified_output1: Option<BufWriter<Box<dyn Write + Send>>>,
    unclassified_output2: Option<BufWriter<Box<dyn Write + Send>>>,
}
impl OutputStreamData {
    fn new() -> Self {
        Self {
            initialized: false,
            kraken_output: Some(BufWriter::new(Box::new(io::stdout()))),
            classified_output1: None,
            classified_output2: None,
            unclassified_output1: None,
            unclassified_output2: None,
        }
    }

    fn initialize(&mut self, opts: &Options, format: SequenceFormat) {
        if self.initialized {
            return;
        }

        // Initialize Classified Outputs
        if let Some(ref filename) = opts.classified_output_filename {
            if opts.paired_end_processing {
                let fields: Vec<&str> = filename.split('#').collect();
                if fields.len() != 2 {
                    eprintln!(
                        "Paired filename format must contain exactly one '#' character: {}",
                        filename
                    );
                    process::exit(1);
                }
                self.classified_output1 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_1{}", fields[0], fields[1])).unwrap(),
                )));
                self.classified_output2 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_2{}", fields[0], fields[1])).unwrap(),
                )));
            } else {
                self.classified_output1 =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            }
        }

        // Initialize Unclassified Outputs
        if let Some(ref filename) = opts.unclassified_output_filename {
            if opts.paired_end_processing {
                let fields: Vec<&str> = filename.split('#').collect();
                if fields.len() != 2 {
                    eprintln!(
                        "Paired filename format must contain exactly one '#' character: {}",
                        filename
                    );
                    process::exit(1);
                }
                self.unclassified_output1 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_1{}", fields[0], fields[1])).unwrap(),
                )));
                self.unclassified_output2 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_2{}", fields[0], fields[1])).unwrap(),
                )));
            } else {
                self.unclassified_output1 =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            }
        }

        // Initialize Kraken Output
        if let Some(ref filename) = opts.kraken_output_filename {
            if filename != "-" {
                self.kraken_output =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            } else {
                self.kraken_output = None;
            }
        }

        self.initialized = true;
    }

    fn write_outputs(&mut self, out_data: &OutputData) {
        if let Some(ref mut kraken_output) = self.kraken_output {
            kraken_output
                .write_all(out_data.kraken_str.as_bytes())
                .unwrap();
        }
        if let Some(ref mut classified_output1) = self.classified_output1 {
            classified_output1
                .write_all(out_data.classified_out1_str.as_bytes())
                .unwrap();
        }
        if let Some(ref mut classified_output2) = self.classified_output2 {
            classified_output2
                .write_all(out_data.classified_out2_str.as_bytes())
                .unwrap();
        }
        if let Some(ref mut unclassified_output1) = self.unclassified_output1 {
            unclassified_output1
                .write_all(out_data.unclassified_out1_str.as_bytes())
                .unwrap();
        }
        if let Some(ref mut unclassified_output2) = self.unclassified_output2 {
            unclassified_output2
                .write_all(out_data.unclassified_out2_str.as_bytes())
                .unwrap();
        }
    }

    fn flush(&mut self) {
        if let Some(ref mut kraken_output) = self.kraken_output {
            kraken_output.flush().unwrap();
        }
        if let Some(ref mut classified_output1) = self.classified_output1 {
            classified_output1.flush().unwrap();
        }
        if let Some(ref mut classified_output2) = self.classified_output2 {
            classified_output2.flush().unwrap();
        }
        if let Some(ref mut unclassified_output1) = self.unclassified_output1 {
            unclassified_output1.flush().unwrap();
        }
        if let Some(ref mut unclassified_output2) = self.unclassified_output2 {
            unclassified_output2.flush().unwrap();
        }
    }
}

#[derive(Default, Clone)]
struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}

impl PartialEq for OutputData {
    fn eq(&self, other: &Self) -> bool {
        self.block_id == other.block_id
    }
}

impl Eq for OutputData {}

impl PartialOrd for OutputData {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OutputData {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.block_id.cmp(&self.block_id)
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut opts = Options::default();
    parse_command_line(&mut opts)?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.num_threads)
        .build_global()?;

    eprintln!("Loading database information...");

    let idx_opts = IndexOptions::load(&opts.options_filename)?;
    opts.use_translated_search = !idx_opts.dna_db;

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping)?;
    let hash = CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping)?;

    eprintln!(" done.");

    let mut stats = ClassificationStats::default();

    let mut outputs = OutputStreamData::default();
    outputs.kraken_output =
        if opts.kraken_output_filename.is_empty() || opts.kraken_output_filename == "-" {
            Some(BufWriter::new(Box::new(io::stdout())))
        } else {
            Some(BufWriter::new(Box::new(File::create(
                &opts.kraken_output_filename,
            )?)))
        };

    let start_time = Instant::now();

    let input_files: Vec<String> = env::args().skip(1).collect();

    let mut taxon_counters: TaxonCounters = HashMap::new();

    if input_files.is_empty() {
        if opts.paired_end_processing && !opts.single_file_pairs {
            eprintln!("Paired-end processing used with no files specified");
            process::exit(1);
        }
        process_files(
            None,
            None,
            &hash,
            &taxonomy,
            &idx_opts,
            &opts,
            &mut stats,
            &mut outputs,
            &mut taxon_counters,
        )?;
    } else {
        let mut i = 0;
        while i < input_files.len() {
            if opts.paired_end_processing && !opts.single_file_pairs {
                if i + 1 == input_files.len() {
                    eprintln!("Paired-end processing used with unpaired file");
                    process::exit(1);
                }
                process_files(
                    Some(&input_files[i]),
                    Some(&input_files[i + 1]),
                    &hash,
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
                    &hash,
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

    let elapsed = start_time.elapsed();
    report_stats(elapsed, &stats);

    if !opts.report_filename.is_empty() {
        if opts.mpa_style_report {
            report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &taxon_counters,
            )?;
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
            )?;
        }
    }

    Ok(())
}

fn report_stats(elapsed: Duration, stats: &ClassificationStats) {
    let seconds = elapsed.as_secs_f64();
    let total_unclassified = stats.total_sequences - stats.total_classified;

    eprint!("\r");
    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1e6,
        seconds,
        stats.total_sequences as f64 / 1e3 / (seconds / 60.0),
        stats.total_bases as f64 / 1e6 / (seconds / 60.0)
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
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &Options,
    stats: &mut ClassificationStats,
    outputs: &mut OutputStreamData,
    total_taxon_counters: &mut TaxonCounters,
) -> Result<(), Box<dyn Error>> {
    let fptr1: Box<dyn BufRead + Send> = match filename1 {
        Some(filename) => Box::new(BufReader::new(File::open(filename)?)),
        None => Box::new(BufReader::new(io::stdin())),
    };
    let fptr2: Option<Box<dyn BufRead + Send>> =
        if opts.paired_end_processing && !opts.single_file_pairs {
            match filename2 {
                Some(filename) => Some(Box::new(BufReader::new(File::open(filename)?))),
                None => {
                    eprintln!("Paired-end processing requires two input files");
                    process::exit(1);
                }
            }
        } else {
            None
        };

    let output_queue = Arc::new(Mutex::new(BinaryHeap::new()));
    let output_lock = Arc::new(Mutex::new(()));

    let mut next_input_block_id = 0u64;
    let next_output_block_id = Arc::new(Mutex::new(0u64));

    let hash = Arc::new(hash.clone());
    let taxonomy = Arc::new(taxonomy.clone());
    let idx_opts = Arc::new(idx_opts.clone());
    let opts = Arc::new(opts.clone());
    let outputs_arc = Arc::new(Mutex::new(outputs.clone()));
    let stats_arc = Arc::new(Mutex::new(stats.clone()));
    let total_taxon_counters_arc = Arc::new(Mutex::new(total_taxon_counters.clone()));

    // Simulating parallel processing with Rayon
    let sequences1 = read_sequences(fptr1)?;
    let sequences2 = if let Some(fptr2) = fptr2 {
        Some(read_sequences(fptr2)?)
    } else {
        None
    };

    let paired = opts.paired_end_processing;

    let chunks: Vec<Vec<Sequence>> = if sequences1.len() > NUM_FRAGMENTS_PER_THREAD {
        sequences1
            .chunks(NUM_FRAGMENTS_PER_THREAD)
            .map(|chunk| chunk.to_vec())
            .collect()
    } else {
        vec![sequences1]
    };

    chunks.into_par_iter().for_each(|chunk| {
        let mut scanner = MinimizerScanner::new(
            idx_opts.k,
            idx_opts.l,
            idx_opts.spaced_seed_mask,
            idx_opts.dna_db,
            idx_opts.toggle_mask,
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
        let mut curr_taxon_counts = TaxonCounters::new();

        let block_id = {
            let id = next_input_block_id;
            next_input_block_id += 1;
            id
        };

        for (i, seq1) in chunk.iter().enumerate() {
            let seq2 = if paired {
                // Get corresponding sequence from sequences2
                sequences2.as_ref().and_then(|s| s.get(i))
            } else {
                None
            };

            if opts.minimum_quality_score > 0 {
                // Mask low-quality bases
                mask_low_quality_bases(seq1, opts.minimum_quality_score);
                if let Some(seq2) = seq2 {
                    mask_low_quality_bases(seq2, opts.minimum_quality_score);
                }
            }

            let call = classify_sequence(
                seq1,
                seq2,
                &mut kraken_oss,
                &hash,
                &taxonomy,
                &idx_opts,
                &opts,
                &mut thread_stats,
                &mut scanner,
                &mut taxa,
                &mut hit_counts,
                &mut translated_frames,
                &mut curr_taxon_counts,
            );

            if call != 0 {
                // Classified output
                let buffer = format!(" kraken:taxid|{}", taxonomy.nodes[&call].external_id);
                let mut seq1_cloned = seq1.clone();
                seq1_cloned.header.push_str(&buffer);
                c1_oss.push_str(&seq1_cloned.to_string());

                if let Some(seq2) = seq2 {
                    let mut seq2_cloned = seq2.clone();
                    seq2_cloned.header.push_str(&buffer);
                    c2_oss.push_str(&seq2_cloned.to_string());
                }
            } else {
                // Unclassified output
                u1_oss.push_str(&seq1.to_string());
                if let Some(seq2) = seq2 {
                    u2_oss.push_str(&seq2.to_string());
                }
            }

            thread_stats.total_sequences += 1;
            thread_stats.total_bases += seq1.seq.len() as u64;
            if let Some(seq2) = seq2 {
                thread_stats.total_bases += seq2.seq.len() as u64;
            }
        }

        // Update global stats
        {
            let mut stats_lock = stats_arc.lock().unwrap();
            stats_lock.total_sequences += thread_stats.total_sequences;
            stats_lock.total_bases += thread_stats.total_bases;
            stats_lock.total_classified += thread_stats.total_classified;
        }

        // Enqueue output data
        let out_data = OutputData {
            block_id,
            kraken_str: kraken_oss,
            classified_out1_str: c1_oss,
            classified_out2_str: c2_oss,
            unclassified_out1_str: u1_oss,
            unclassified_out2_str: u2_oss,
        };

        {
            let mut queue = output_queue.lock().unwrap();
            queue.push(Reverse(out_data));
        }

        // Update taxon counters
        {
            let mut total_counters = total_taxon_counters_arc.lock().unwrap();
            for (tax_id, counter) in curr_taxon_counts {
                total_counters
                    .entry(tax_id)
                    .or_insert_with(TaxonCount::default)
                    .merge(&counter);
            }
        }

        // Output loop
        loop {
            let should_output = {
                let queue = output_queue.lock().unwrap();
                if let Some(Reverse(front)) = queue.peek() {
                    front.block_id == *next_output_block_id.lock().unwrap()
                } else {
                    false
                }
            };

            if should_output {
                let out_data = {
                    let mut queue = output_queue.lock().unwrap();
                    queue.pop().unwrap().0
                };

                let _lock = output_lock.lock().unwrap();

                // Output data
                let mut outputs_lock = outputs_arc.lock().unwrap();
                if !outputs_lock.initialized {
                    initialize_outputs(&opts, &mut outputs_lock, seq1.format);
                }
                if let Some(ref mut kraken_output) = outputs_lock.kraken_output {
                    kraken_output
                        .write_all(out_data.kraken_str.as_bytes())
                        .unwrap();
                }
                if let Some(ref mut classified_output1) = outputs_lock.classified_output1 {
                    classified_output1
                        .write_all(out_data.classified_out1_str.as_bytes())
                        .unwrap();
                }
                if let Some(ref mut classified_output2) = outputs_lock.classified_output2 {
                    classified_output2
                        .write_all(out_data.classified_out2_str.as_bytes())
                        .unwrap();
                }
                if let Some(ref mut unclassified_output1) = outputs_lock.unclassified_output1 {
                    unclassified_output1
                        .write_all(out_data.unclassified_out1_str.as_bytes())
                        .unwrap();
                }
                if let Some(ref mut unclassified_output2) = outputs_lock.unclassified_output2 {
                    unclassified_output2
                        .write_all(out_data.unclassified_out2_str.as_bytes())
                        .unwrap();
                }

                let mut next_output = next_output_block_id.lock().unwrap();
                *next_output += 1;
            } else {
                break;
            }
        }
    });

    // Flush outputs
    if let Some(ref mut kraken_output) = outputs.kraken_output {
        kraken_output.flush()?;
    }
    if let Some(ref mut classified_output1) = outputs.classified_output1 {
        classified_output1.flush()?;
    }
    if let Some(ref mut classified_output2) = outputs.classified_output2 {
        classified_output2.flush()?;
    }
    if let Some(ref mut unclassified_output1) = outputs.unclassified_output1 {
        unclassified_output1.flush()?;
    }
    if let Some(ref mut unclassified_output2) = outputs.unclassified_output2 {
        unclassified_output2.flush()?;
    }

    Ok(())
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
    hit_counts: &mut TaxonCounts,
    tx_frames: &mut Vec<String>,
    curr_taxon_counts: &mut TaxonCounters,
) -> TaxId {
    let mut call = 0;
    taxa.clear();
    hit_counts.clear();
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups = 0;

    'outer: for mate_num in 0..2 {
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        if opts.use_translated_search {
            let seq = if mate_num == 0 { &dna.seq } else { &dna2.seq };
            translate_to_all_frames(seq, tx_frames);
        }

        for frame_idx in 0..frame_ct {
            if opts.use_translated_search {
                scanner.load_sequence(&tx_frames[frame_idx]);
            } else {
                let seq = if mate_num == 0 { &dna.seq } else { &dna2.seq };
                scanner.load_sequence(seq);
            }

            let mut last_minimizer = UINT64_MAX;
            let mut last_taxon = MATE_PAIR_BORDER_TAXON; // Initialize to a value

            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon = if scanner.is_ambiguous() {
                    AMBIGUOUS_SPAN_TAXON
                } else {
                    if minimizer != last_minimizer {
                        let mut skip_lookup = false;
                        if idx_opts.minimum_acceptable_hash_value > 0 {
                            if murmur_hash3(minimizer) < idx_opts.minimum_acceptable_hash_value {
                                skip_lookup = true;
                            }
                        }

                        let taxon = if !skip_lookup {
                            hash.get(minimizer).unwrap_or(0)
                        } else {
                            0
                        };

                        last_taxon = taxon;
                        last_minimizer = minimizer;

                        if taxon != 0 {
                            minimizer_hit_groups += 1;
                            curr_taxon_counts
                                .entry(taxon)
                                .or_insert_with(TaxonCount::default)
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
                        break 'outer; // Break out of all loops
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
    }

    let mut total_kmers = taxa.len();
    if opts.paired_end_processing {
        total_kmers -= 1; // account for the mate pair marker
    }
    if opts.use_translated_search {
        total_kmers -= if opts.paired_end_processing { 4 } else { 2 }; // account for reading frame markers
    }

    call = resolve_tree(hit_counts, taxonomy, total_kmers, opts);

    // Void a call made by too few minimizer groups
    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    if call != 0 {
        stats.total_classified += 1;
        curr_taxon_counts
            .entry(call)
            .or_insert_with(TaxonCount::default)
            .increment_read_count();
    }

    if call != 0 {
        koss.push_str("C\t");
    } else {
        koss.push_str("U\t");
    }

    if !opts.paired_end_processing {
        koss.push_str(&format!("{}\t", dna.id));
    } else {
        let trimmed_id = trim_pair_info(&dna.id);
        koss.push_str(&format!("{}\t", trimmed_id));
    }

    let ext_call = taxonomy.nodes()[call as usize].external_id;
    if opts.print_scientific_name {
        let name = if call != 0 {
            let node = &taxonomy.nodes()[call as usize];
            let name_offset = node.name_offset;
            Some(&taxonomy.name_data()[name_offset..])
        } else {
            None
        };
        let name_str = name.unwrap_or("unclassified");
        koss.push_str(&format!("{} (taxid {})", name_str, ext_call));
    } else {
        koss.push_str(&format!("{}", ext_call));
    }

    koss.push('\t');
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
            add_hitlist_string(koss, taxa, taxonomy);
        }
    }

    koss.push('\n');

    call
}

fn resolve_tree(
    hit_counts: &HashMap<TaxId, u32>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (opts.confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for &taxon in hit_counts.keys() {
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
    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);
    // We probably have a call w/o required support (unless LCA resolved tie)
    while max_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts.iter() {
            if taxonomy.is_a_ancestor_of_b(max_taxon, taxon) {
                max_score += count;
            }
        }
        // Score is now sum of hits at max_taxon and within max_taxon clade
        if max_score >= required_score {
            return max_taxon;
        } else {
            // Run up tree until confidence threshold is met
            // Run off tree if required score isn't met
            max_taxon = taxonomy.nodes[&max_taxon].parent_id;
        }
    }

    max_taxon
}

fn add_hitlist_string(koss: &mut String, taxa: &[TaxId], taxonomy: &Taxonomy) {
    if taxa.is_empty() {
        return;
    }

    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in &taxa[1..] {
        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    koss.push_str(&format!("A:{} ", code_count));
                } else {
                    if let Some(node) = taxonomy.nodes.get(&last_code) {
                        let ext_code = node.external_id;
                        koss.push_str(&format!("{}:{} ", ext_code, code_count));
                    }
                }
            } else {
                // Mate pair/reading frame marker
                let marker = if last_code == MATE_PAIR_BORDER_TAXON {
                    "|:| "
                } else {
                    "-:- "
                };
                koss.push_str(marker);
            }
            code_count = 1;
            last_code = code;
        }
    }

    // Handle the last code
    if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
        if last_code == AMBIGUOUS_SPAN_TAXON {
            koss.push_str(&format!("A:{}", code_count));
        } else {
            if let Some(node) = taxonomy.nodes.get(&last_code) {
                let ext_code = node.external_id;
                koss.push_str(&format!("{}:{}", ext_code, code_count));
            }
        }
    } else {
        let marker = if last_code == MATE_PAIR_BORDER_TAXON {
            "|:|"
        } else {
            "-:-"
        };
        koss.push_str(marker);
    }
}

fn initialize_outputs(opts: &Options, outputs: &mut OutputStreamData, format: SequenceFormat) {
    if !outputs.initialized {
        // Classified Outputs
        if let Some(ref filename) = opts.classified_output_filename {
            if opts.paired_end_processing {
                let fields: Vec<&str> = filename.split('#').collect();
                if fields.len() != 2 {
                    eprintln!(
                        "Paired filename format must contain exactly one '#' character: {}",
                        filename
                    );
                    process::exit(1);
                }
                outputs.classified_output1 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_1{}", fields[0], fields[1])).unwrap(),
                )));
                outputs.classified_output2 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_2{}", fields[0], fields[1])).unwrap(),
                )));
            } else {
                outputs.classified_output1 =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            }
        }

        // Unclassified Outputs
        if let Some(ref filename) = opts.unclassified_output_filename {
            if opts.paired_end_processing {
                let fields: Vec<&str> = filename.split('#').collect();
                if fields.len() != 2 {
                    eprintln!(
                        "Paired filename format must contain exactly one '#' character: {}",
                        filename
                    );
                    process::exit(1);
                }
                outputs.unclassified_output1 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_1{}", fields[0], fields[1])).unwrap(),
                )));
                outputs.unclassified_output2 = Some(BufWriter::new(Box::new(
                    File::create(format!("{}_2{}", fields[0], fields[1])).unwrap(),
                )));
            } else {
                outputs.unclassified_output1 =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            }
        }

        // Kraken Output
        if let Some(ref filename) = opts.kraken_output_filename {
            if filename != "-" {
                outputs.kraken_output =
                    Some(BufWriter::new(Box::new(File::create(filename).unwrap())));
            } else {
                outputs.kraken_output = None;
            }
        }

        outputs.initialized = true;
    }
}

fn mask_low_quality_bases(dna: &mut Sequence, minimum_quality_score: usize) {
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
        process::exit(1);
    }
    let quals = dna.quals.as_bytes();
    let mut seq = dna.seq.clone().into_bytes();
    for (i, &q) in quals.iter().enumerate() {
        if q - b'!' < minimum_quality_score as u8 {
            seq[i] = b'x';
        }
    }
    dna.seq = String::from_utf8(seq).unwrap();
}

fn trim_pair_info(id: &str) -> String {
    if id.len() <= 2 {
        return id.to_string();
    }
    if id.ends_with("/1") || id.ends_with("/2") {
        id[..id.len() - 2].to_string()
    } else {
        id.to_string()
    }
}

fn parse_command_line(opts: &mut Options) -> Result<(), Box<dyn Error>> {
    let matches = Command::new("classify")
        .version("1.0")
        .about("Kraken 2 taxonomic sequence classification system")
        .arg(
            Arg::new("index_filename")
                .short('H')
                .long("index")
                .action(ArgAction::Set)
                .required(true)
                .help("Kraken 2 index filename"),
        )
        .arg(
            Arg::new("taxonomy_filename")
                .short('t')
                .long("taxonomy")
                .action(ArgAction::Set)
                .required(true)
                .help("Kraken 2 taxonomy filename"),
        )
        .arg(
            Arg::new("options_filename")
                .short('o')
                .long("options")
                .action(ArgAction::Set)
                .required(true)
                .help("Kraken 2 options filename"),
        )
        .arg(
            Arg::new("quick_mode")
                .short('q')
                .long("quick")
                .action(ArgAction::SetTrue)
                .help("Quick mode"),
        )
        .arg(
            Arg::new("use_memory_mapping")
                .short('M')
                .long("memory-map")
                .action(ArgAction::SetTrue)
                .help("Use memory mapping to access hash & taxonomy"),
        )
        .arg(
            Arg::new("confidence_threshold")
                .short('T')
                .long("confidence")
                .action(ArgAction::Set)
                .help("Confidence score threshold (default 0)"),
        )
        .arg(
            Arg::new("num_threads")
                .short('p')
                .long("threads")
                .action(ArgAction::Set)
                .help("Number of threads (default 1)"),
        )
        .arg(
            Arg::new("minimum_quality_score")
                .short('Q')
                .long("min-quality")
                .action(ArgAction::Set)
                .help("Minimum quality score (FASTQ only, default 0)"),
        )
        .arg(
            Arg::new("paired_end_processing")
                .short('P')
                .long("paired")
                .action(ArgAction::SetTrue)
                .help("Process pairs of reads"),
        )
        .arg(
            Arg::new("single_file_pairs")
                .short('S')
                .long("single-file-pairs")
                .action(ArgAction::SetTrue)
                .help("Process pairs with mates in same file"),
        )
        .arg(
            Arg::new("report_filename")
                .short('R')
                .long("report")
                .action(ArgAction::Set)
                .help("Print report to filename"),
        )
        .arg(
            Arg::new("mpa_style_report")
                .short('m')
                .long("mpa-report")
                .action(ArgAction::SetTrue)
                .help("In combination with -R, use mpa-style report"),
        )
        .arg(
            Arg::new("report_zero_counts")
                .short('z')
                .long("report-zero")
                .action(ArgAction::SetTrue)
                .help("In combination with -R, report taxa with zero count"),
        )
        .arg(
            Arg::new("print_scientific_name")
                .short('n')
                .long("scientific-name")
                .action(ArgAction::SetTrue)
                .help("Print scientific name instead of taxid in Kraken output"),
        )
        .arg(
            Arg::new("minimum_hit_groups")
                .short('g')
                .long("min-hit-groups")
                .action(ArgAction::Set)
                .help("Minimum number of hit groups needed for call"),
        )
        .arg(
            Arg::new("classified_output_filename")
                .short('C')
                .long("classified")
                .action(ArgAction::Set)
                .help("Filename/format to have classified sequences"),
        )
        .arg(
            Arg::new("unclassified_output_filename")
                .short('U')
                .long("unclassified")
                .action(ArgAction::Set)
                .help("Filename/format to have unclassified sequences"),
        )
        .arg(
            Arg::new("kraken_output_filename")
                .short('O')
                .long("kraken-output")
                .action(ArgAction::Set)
                .help("Output file for normal Kraken output"),
        )
        .arg(
            Arg::new("report_kmer_data")
                .short('K')
                .long("report-kmer")
                .action(ArgAction::SetTrue)
                .help("In combination with -R, provide minimizer information in report"),
        )
        .get_matches();

    opts.index_filename = matches
        .get_one::<String>("index_filename")
        .unwrap()
        .to_string();
    opts.taxonomy_filename = matches
        .get_one::<String>("taxonomy_filename")
        .unwrap()
        .to_string();
    opts.options_filename = matches
        .get_one::<String>("options_filename")
        .unwrap()
        .to_string();

    opts.quick_mode = matches.get_flag("quick_mode");
    opts.use_memory_mapping = matches.get_flag("use_memory_mapping");
    opts.paired_end_processing = matches.get_flag("paired_end_processing");
    opts.single_file_pairs = matches.get_flag("single_file_pairs");
    opts.mpa_style_report = matches.get_flag("mpa_style_report");
    opts.report_zero_counts = matches.get_flag("report_zero_counts");
    opts.print_scientific_name = matches.get_flag("print_scientific_name");
    opts.report_kmer_data = matches.get_flag("report_kmer_data");

    if let Some(value) = matches.get_one::<String>("confidence_threshold") {
        opts.confidence_threshold = value.parse()?;
        if opts.confidence_threshold < 0.0 || opts.confidence_threshold > 1.0 {
            eprintln!("Confidence threshold must be in [0, 1]");
            process::exit(1);
        }
    }

    if let Some(value) = matches.get_one::<String>("num_threads") {
        opts.num_threads = value.parse()?;
        if opts.num_threads == 0 {
            eprintln!("Number of threads can't be less than 1");
            process::exit(1);
        }
    }

    if let Some(value) = matches.get_one::<String>("minimum_quality_score") {
        opts.minimum_quality_score = value.parse()?;
    }

    if let Some(value) = matches.get_one::<String>("minimum_hit_groups") {
        opts.minimum_hit_groups = value.parse()?;
    }

    if let Some(value) = matches.get_one::<String>("report_filename") {
        opts.report_filename = Some(value.to_string());
    }

    if let Some(value) = matches.get_one::<String>("classified_output_filename") {
        opts.classified_output_filename = Some(value.to_string());
    }

    if let Some(value) = matches.get_one::<String>("unclassified_output_filename") {
        opts.unclassified_output_filename = Some(value.to_string());
    }

    if let Some(value) = matches.get_one::<String>("kraken_output_filename") {
        opts.kraken_output_filename = Some(value.to_string());
    }

    if opts.mpa_style_report && opts.report_filename.is_none() {
        eprintln!("-m requires -R to be used");
        process::exit(1);
    }

    Ok(())
}

fn read_sequences<R: BufRead + Send>(reader: R) -> Result<Vec<Sequence>, Box<dyn Error>> {
    let mut sequences = Vec::new();
    let mut lines = reader.lines();
    while let Some(line) = lines.next() {
        let id = line?;
        if let Some(seq_line) = lines.next() {
            let seq = seq_line?;
            sequences.push(Sequence {
                id,
                seq,
                header: String::new(),
                quals: String::new(),
                format: SequenceFormat::Fasta, // Adjust based on actual format
            });
        }
    }
    Ok(sequences)
}

fn translate_to_all_frames(seq: &str, tx_frames: &mut Vec<String>) {
    // Implement translation to all six reading frames
    // For simplicity, we'll push the original sequence multiple times
    // Replace with actual translation logic
    tx_frames.clear();
    for frame in 0..6 {
        tx_frames.push(format!("Frame {}: {}", frame + 1, seq));
    }
}

// ... Implement other necessary structs and functions (e.g., Sequence, Taxonomy, MinimizerScanner, etc.)
