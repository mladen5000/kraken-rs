use crossbeam_channel::bounded;
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use clap::{Arg, Command};
use lazy_static::lazy_static;

// Assuming these are defined elsewhere
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{IndexOptions, TaxId};
use crate::kv_store::murmur_hash3;
use crate::mmscanner::MinimizerScanner;
use crate::readcounts::{self, HyperLogLogPlusMinus, ReadCounts};
use crate::reports::{report_kraken_style, report_mpa_style};
use crate::seqreader::{Sequence, SequenceFormat};
use crate::taxonomy::Taxonomy;

lazy_static! {
    static ref GLOBAL_TAXON_COUNTERS: Mutex<HashMap<TaxId, ReadCounts<readcounts::HyperLogLogPlusMinus>>> =
        Mutex::new(HashMap::new());
}

const NUM_FRAGMENTS_PER_THREAD: usize = 10000;
const MATE_PAIR_BORDER_TAXON: TaxId = TaxId::MAX;
const READING_FRAME_BORDER_TAXON: TaxId = TaxId::MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxId = TaxId::MAX - 2;

#[derive(Clone)]
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

#[derive(Default)]
struct OutputStreamData {
    initialized: bool,
    printing_sequences: bool,
    classified_output1: Option<Mutex<Box<dyn io::Write + Send>>>,
    classified_output2: Option<Mutex<Box<dyn io::Write + Send>>>,
    unclassified_output1: Option<Mutex<Box<dyn io::Write + Send>>>,
    unclassified_output2: Option<Mutex<Box<dyn io::Write + Send>>>,
    kraken_output: Option<Mutex<Box<dyn io::Write + Send>>>,
}

#[derive(Debug, Clone)]
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

impl PartialEq for OutputData {
    fn eq(&self, other: &Self) -> bool {
        self.block_id == other.block_id
    }
}

impl Eq for OutputData {}

fn main() -> io::Result<()> {
    let mut opts = parse_command_line();

    println!("Loading database information...");
    let idx_opts = load_index_options(&opts.options_filename)?;
    opts.use_translated_search = !idx_opts.dna_db;

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename, opts.use_memory_mapping)?;
    // 10000 entries, 24, 8
    // let hash = CompactHashTable::new(&opts.index_filename, opts.use_memory_mapping, &idx_opts)?;
    let hash = CompactHashTable::new(10000, 24, 8);

    println!("done.");

    let stats = Arc::new(Mutex::new(ClassificationStats::default()));
    let outputs = Arc::new(Mutex::new(OutputStreamData {
        initialized: false,
        printing_sequences: false,
        classified_output1: None,
        classified_output2: None,
        unclassified_output1: None,
        unclassified_output2: None,
        kraken_output: Some(Mutex::new(Box::new(io::stdout()))),
    }));

    let start_time = Instant::now();

    // Extract input filenames from command-line arguments
    let args: Vec<String> = std::env::args().collect();
    let filename1 = if opts.single_file_pairs {
        opts.taxonomy_filename.as_str()
    } else if args.len() > 1 {
        args[1].as_str()
    } else {
        ""
    };

    let filename2 = if opts.paired_end_processing && args.len() > 2 {
        if args.len() > 2 {
            Some(args[2].as_str())
        } else {
            None
        }
    } else {
        None
    };

    if args.len() <= 1 {
        if opts.paired_end_processing && !opts.single_file_pairs {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "paired end processing used with no files specified",
            ));
        }
        process_files(
            Some(filename1),
            filename2,
            Arc::new(hash?),
            opts.taxonomy_filename.clone(),
            Arc::new(idx_opts),
            Arc::new(opts.clone()),
            Arc::new(Mutex::new(ClassificationStats::default())),
            Arc::new(Mutex::new(OutputStreamData::default())),
        )?;
    } else {
        process_files(
            Some(filename1),
            filename2,
            Arc::new(hash?),
            opts.taxonomy_filename.clone(),
            Arc::new(idx_opts),
            Arc::new(opts.clone()),
            Arc::clone(&stats),
            Arc::clone(&outputs),
        )?;
    }

    let elapsed = start_time.elapsed();

    // Lock the Mutex to access ClassificationStats
    let stats_guard = stats.lock().unwrap();

    // Call report_stats with accessed fields
    report_stats(elapsed, &*stats_guard);

    // Initialize taxon counters
    let taxon_counters = GLOBAL_TAXON_COUNTERS.lock().unwrap();
    if !opts.report_filename.is_empty() {
        // Lock the Mutex to access ClassificationStats
        let stats_guard = stats.lock().unwrap();

        // Calculate total_unclassified using the locked data
        let total_unclassified = stats_guard.total_sequences - stats_guard.total_classified;

        if opts.mpa_style_report {
            report_mpa_style(
                &opts.report_filename,
                opts.report_zero_counts,
                &taxonomy,
                &*taxon_counters,
            )?;
        } else {
            // Lock the Mutex to access ClassificationStats
            let stats_guard = stats.lock().unwrap();

            // Calculate total_unclassified using the locked data
            let total_unclassified = stats_guard.total_sequences - stats_guard.total_classified;

            // Call report_kraken_style with accessed fields
            report_kraken_style(
                &opts.report_filename,
                opts.report_zero_counts,
                opts.report_kmer_data,
                &taxonomy,
                &*taxon_counters,
                stats_guard.total_sequences,
                total_unclassified,
            )?;
        }
    }

    Ok(())
}

pub fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: Arc<CompactHashTable>,
    taxonomy_filename: String,
    idx_opts: Arc<IndexOptions>,
    opts: Arc<Options>,
    stats: Arc<Mutex<ClassificationStats>>,
    outputs: Arc<Mutex<OutputStreamData>>,
) -> io::Result<()> {
    let (sender, receiver) = bounded(opts.num_threads * 2);
    let reader1: Box<dyn BufRead> = match filename1 {
        Some(f) => Box::new(BufReader::new(File::open(f)?)),
        None => Box::new(io::stdin().lock()),
    };
    let reader2 = filename2.map(|f| BufReader::new(File::open(f).unwrap()));

    let next_input_block_id = Arc::new(Mutex::new(0u64));
    let next_output_block_id = Arc::new(Mutex::new(0u64));
    let output_queue = Arc::new(Mutex::new(BinaryHeap::new()));

    let reader_opts = Arc::clone(&opts);
    let reader_thread = thread::spawn(move || -> io::Result<()> {
        loop {
            let sequences1 = Vec::new();
            let sequences2 = if reader_opts.paired_end_processing {
                Some(Vec::new())
            } else {
                None
            };
            // Read sequences (similar to your implementation)
            // ...

            if sequences1.is_empty() {
                break;
            }

            let block_id = {
                let mut id = next_input_block_id.lock().unwrap();
                let current_id = *id;
                *id += 1;
                current_id
            };

            sender.send((block_id, sequences1, sequences2)).unwrap();
        }
        Ok(())
    });

    let worker_threads: Vec<_> = (0..opts.num_threads)
        .map(|_| {
            let receiver = receiver.clone();
            let hash = Arc::clone(&hash);
            let taxonomy_filename = taxonomy_filename.clone();
            let opts = Arc::clone(&opts);
            let idx_opts = Arc::clone(&idx_opts);
            let stats = Arc::clone(&stats);
            let outputs = Arc::clone(&outputs);
            let next_output_block_id = Arc::clone(&next_output_block_id);
            let output_queue = Arc::clone(&output_queue);

            thread::spawn(move || {
                let tax = Taxonomy::new(&taxonomy_filename, opts.use_memory_mapping).unwrap();
                let mut thread_stats = ClassificationStats::default();
                let mut thread_taxon_counters: HashMap<
                    TaxId,
                    ReadCounts<readcounts::HyperLogLogPlusMinus>,
                > = HashMap::new();

                let mut scanner = MinimizerScanner::new(
                    idx_opts.k as isize,
                    idx_opts.l as isize,
                    idx_opts.spaced_seed_mask,
                    idx_opts.dna_db,
                    idx_opts.toggle_mask,
                    idx_opts.revcom_version as i32,
                );

                while let Ok((block_id, sequences1, sequences2)) = receiver.recv() {
                    let kraken_oss = String::new();
                    let c1_oss = String::new();
                    let c2_oss = String::new();
                    let u1_oss = String::new();
                    let u2_oss = String::new();

                    for (i, seq1) in sequences1.iter().enumerate() {
                        let seq2 = sequences2.as_ref().and_then(|seqs| seqs.get(i));

                        let call = classify_sequence(
                            seq1,
                            seq2,
                            &hash,
                            &tax,
                            &idx_opts,
                            &opts,
                            &mut thread_stats,
                            &mut scanner,
                        );

                        thread_stats.total_sequences += 1;
                        thread_stats.total_bases += seq1.seq.len() as u64;
                        if let Some(seq2) = seq2 {
                            thread_stats.total_bases += seq2.seq.len() as u64;
                        }
                        thread_taxon_counters
                            .entry(call)
                            .or_insert_with(ReadCounts::default)
                            .increment_read_count();
                    }

                    let out_data = OutputData {
                        block_id,
                        kraken_str: kraken_oss,
                        classified_out1_str: c1_oss,
                        classified_out2_str: c2_oss,
                        unclassified_out1_str: u1_oss,
                        unclassified_out2_str: u2_oss,
                    };

                    // Add output data to queue
                    let mut queue = output_queue.lock().unwrap();
                    queue.push(out_data);

                    // Process output queue
                    while let Some(out_data) = queue.peek() {
                        let mut next_output_block_id = next_output_block_id.lock().unwrap();
                        if out_data.block_id == *next_output_block_id {
                            let out_data = queue.pop().unwrap();
                            drop(queue);
                            write_output(&out_data, &outputs).unwrap();
                            *next_output_block_id += 1;
                            queue = output_queue.lock().unwrap();
                        } else {
                            break;
                        }
                    }

                    // Update taxon counters
                    // This would be implemented similar to the C++ version
                    // ...
                }
                // Update global stats
                let mut global_stats = stats.lock().unwrap();
                global_stats.total_sequences += thread_stats.total_sequences;
                global_stats.total_bases += thread_stats.total_bases;
                global_stats.total_classified += thread_stats.total_classified;
                // Update global taxon counters (you'll need to implement this part)
                update_global_taxon_counters(&thread_taxon_counters);
            })
        })
        .collect();

    // Wait for all threads to complete
    reader_thread.join().unwrap()?;
    for thread in worker_threads {
        thread.join().unwrap();
    }

    // Process any remaining items in the output queue
    let mut queue = output_queue.lock().unwrap();
    while let Some(out_data) = queue.pop() {
        write_output(&out_data, &outputs)?;
    }

    Ok(())
}

fn classify_sequence<'a, 'b>(
    dna: &'b Sequence,
    dna2: Option<&'b Sequence>,
    hash: &'b CompactHashTable,
    tax: &'b Taxonomy,
    idx_opts: &'b IndexOptions,
    opts: &'b Options,
    stats: &'b mut ClassificationStats,
    scanner: &'b mut MinimizerScanner<'a>,
) -> u64
where
    'a: 'b,
{
    // Initialize vectors and counters
    let mut taxa = Vec::new();
    let mut hit_counts = HashMap::new();
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups = 0usize;

    // Process each mate (for paired-end sequences)
    for mate_num in 0..2 {
        // Break if not processing paired-end sequences and on second iteration
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        // Select the appropriate sequence based on mate number
        let seq = if mate_num == 0 {
            &dna.seq
        } else {
            match dna2 {
                Some(seq) => &seq.seq,
                None => break, // Exit the loop if dna2 is None
            }
        };

        // Process the sequence
        {
            // Translate sequences if using translated search
            // if opts.use_translated_search {
            //     let tx_frames = translate_to_all_frames(seq, frame_ct);
            //     for frame_idx in 0..frame_ct {
            //         scanner.load_sequence(&tx_frames[frame_idx], frame_idx, frame_ct);
            //     }
            // } else {
            //     scanner.load_sequence(seq, 0, frame_ct);
            // }

            let mut last_minimizer = 0u64;

            // Process each minimizer in the sequence
            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon = if scanner.is_ambiguous() {
                    AMBIGUOUS_SPAN_TAXON
                } else if minimizer != last_minimizer {
                    // Determine taxon based on minimizer
                    let taxon = if idx_opts.minimum_acceptable_hash_value == 0
                        || murmur_hash3(minimizer) >= idx_opts.minimum_acceptable_hash_value
                    {
                        hash.get(minimizer)
                    } else {
                        0
                    };

                    if taxon != 0 {
                        minimizer_hit_groups += 1;
                        // TODO: Update taxon counters here
                    }

                    last_minimizer = minimizer;
                    taxon
                } else {
                    continue;
                };

                if taxon != 0 {
                    // Return early if in quick mode and minimum hit groups reached
                    if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                        return taxon;
                    }
                    *hit_counts.entry(taxon).or_insert(0) += 1;
                }

                taxa.push(taxon);
            }
        }

        // Add frame and mate pair border taxa
        if opts.use_translated_search {
            taxa.push(READING_FRAME_BORDER_TAXON);
        }
        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    // Calculate total k-mers
    let total_kmers = {
        let mut count = taxa.len();
        if opts.paired_end_processing {
            count -= 1;
        }
        if opts.use_translated_search {
            count -= if opts.paired_end_processing { 4 } else { 2 };
        }
        count
    };

    // Resolve the taxonomy tree and make the classification call

    // Generate Kraken output string here (implementation needed)

    {
        let call = resolve_tree(&hit_counts, tax, total_kmers, opts);
        if call == 0 || minimizer_hit_groups < opts.minimum_hit_groups {
            0
        } else {
            stats.total_classified += 1;
            // TODO: Update taxon counters for the call
            call
        }
    }
}

fn resolve_tree(
    hit_counts: &HashMap<TaxId, u32>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    opts: &Options,
) -> TaxId {
    let mut max_taxon = 0;
    let mut max_score = 0u32;
    let required_score = (opts.confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for (&taxon, &count) in hit_counts {
        let score: u32 = hit_counts
            .iter()
            .filter(|(&t, _)| taxonomy.is_a_ancestor_of_b(t, taxon))
            .map(|(_, &c)| c)
            .sum();

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
        max_score = hit_counts
            .iter()
            .filter(|(&t, _)| taxonomy.is_a_ancestor_of_b(max_taxon, t))
            .map(|(_, &c)| c)
            .sum();

        if max_score >= required_score {
            return max_taxon;
        } else {
            max_taxon = taxonomy.get_parent(max_taxon).unwrap_or(0);
        }
    }

    max_taxon
}
fn add_hitlist_string(oss: &mut String, taxa: &[TaxId], taxonomy: &Taxonomy) {
    use std::fmt::Write;
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in taxa.iter().skip(1) {
        if code == last_code {
            code_count += 1;
        } else {
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    write!(oss, "A:{} ", code_count).unwrap();
                } else {
                    let ext_code = taxonomy.nodes[last_code as usize].external_id;
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
            write!(oss, "A:{}", code_count).unwrap();
        } else {
            let ext_code = taxonomy.nodes[last_code as usize].external_id;
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

fn initialize_outputs(
    opts: &Options,
    outputs: &Mutex<OutputStreamData>,
    format: SequenceFormat,
) -> io::Result<()> {
    let mut outputs = outputs.lock().unwrap();
    if !outputs.initialized {
        if !opts.classified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<&str> = opts.classified_output_filename.split('#').collect();
                if fields.len() != 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Paired filename format must contain exactly one '#' character",
                    ));
                }
                outputs.classified_output1 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(format!("{}_1{}", fields[0], fields[1]))?,
                ))));
                outputs.classified_output2 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(format!("{}_2{}", fields[0], fields[1]))?,
                ))));
            } else {
                outputs.classified_output1 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(&opts.classified_output_filename)?,
                ))));
            }
            outputs.printing_sequences = true;
        }

        if !opts.unclassified_output_filename.is_empty() {
            if opts.paired_end_processing {
                let fields: Vec<&str> = opts.unclassified_output_filename.split('#').collect();
                if fields.len() != 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Paired filename format must contain exactly one '#' character",
                    ));
                }
                outputs.unclassified_output1 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(format!("{}_1{}", fields[0], fields[1]))?,
                ))));
                outputs.unclassified_output2 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(format!("{}_1{}", fields[0], fields[1]))?,
                ))));
            } else {
                outputs.unclassified_output1 = Some(Mutex::new(Box::new(BufWriter::new(
                    File::create(&opts.unclassified_output_filename)?,
                ))));
            }
            outputs.printing_sequences = true;
        }

        if !opts.kraken_output_filename.is_empty() {
            if opts.kraken_output_filename != "-" {
                outputs.kraken_output = Some(Mutex::new(Box::new(BufWriter::new(File::create(
                    &opts.kraken_output_filename,
                )?))));
            } else {
                outputs.kraken_output = Some(Mutex::new(Box::new(io::stdout())));
            }
        }

        if !opts.kraken_output_filename.is_empty() && opts.kraken_output_filename != "-" {
            outputs.kraken_output = Some(Mutex::new(Box::new(BufWriter::new(File::create(
                &opts.kraken_output_filename,
            )?))));
        }

        outputs.initialized = true;
    }
    Ok(())
}

fn mask_low_quality_bases(
    dna: &mut Sequence,
    minimum_quality_score: u8,
) -> Result<(), &'static str> {
    if dna.format != SequenceFormat::Fastq {
        return Ok(());
    }

    if dna.seq.len() != dna.quals.len() {
        return Err("Sequence length does not match Quality string length");
    }

    for (base, qual) in unsafe {
        dna.seq
            .as_bytes_mut()
            .iter_mut()
            .zip(dna.quals.as_bytes().iter())
    } {
        if (*qual - b'!') < minimum_quality_score {
            *base = b'x';
        }
    }

    Ok(())
}

fn parse_command_line() -> Options {
    let matches = Command::new("classify")
        .arg(
            Arg::new("index")
                .short('H')
                .long("index")
                .value_name("FILE")
                .help("Kraken 2 index filename")
                .required(true),
        )
        .arg(
            Arg::new("taxonomy")
                .short('t')
                .long("taxonomy")
                .value_name("FILE")
                .help("Kraken 2 taxonomy filename")
                .required(true),
        )
        .arg(
            Arg::new("options")
                .short('o')
                .long("options")
                .value_name("FILE")
                .help("Kraken 2 options filename")
                .required(true),
        )
        .arg(
            Arg::new("quick")
                .short('q')
                .long("quick")
                .help("Quick mode"),
        )
        .arg(
            Arg::new("memory-mapping")
                .short('M')
                .long("memory-mapping")
                .help("Use memory mapping to access hash & taxonomy"),
        )
        .arg(
            Arg::new("confidence")
                .short('T')
                .long("confidence")
                .value_name("NUM")
                .help("Confidence score threshold (def. 0)"),
        )
        .arg(
            Arg::new("threads")
                .short('p')
                .long("threads")
                .value_name("NUM")
                .help("Number of threads (def. 1)"),
        )
        .arg(
            Arg::new("min-quality")
                .short('Q')
                .long("min-quality")
                .value_name("NUM")
                .help("Minimum quality score (FASTQ only, def. 0)"),
        )
        .arg(
            Arg::new("paired")
                .short('P')
                .long("paired")
                .help("Process pairs of reads"),
        )
        .arg(
            Arg::new("single-file-pairs")
                .short('S')
                .long("single-file-pairs")
                .help("Process pairs with mates in same file"),
        )
        .arg(
            Arg::new("report")
                .short('R')
                .long("report")
                .value_name("FILE")
                .help("Print report to filename"),
        )
        .arg(
            Arg::new("mpa-style")
                .short('m')
                .long("mpa-style")
                .help("In comb. w/ -R, use mpa-style report"),
        )
        .arg(
            Arg::new("report-zero-counts")
                .short('z')
                .long("report-zero-counts")
                .help("In comb. w/ -R, report taxa w/ 0 count"),
        )
        .arg(
            Arg::new("scientific-name")
                .short('n')
                .long("scientific-name")
                .help("Print scientific name instead of taxid in Kraken output"),
        )
        .arg(
            Arg::new("hit-groups")
                .short('g')
                .long("hit-groups")
                .value_name("NUM")
                .help("Minimum number of hit groups needed for call"),
        )
        .arg(
            Arg::new("classified-output")
                .short('C')
                .long("classified-output")
                .value_name("FILE")
                .help("Filename/format to have classified sequences"),
        )
        .arg(
            Arg::new("unclassified-output")
                .short('U')
                .long("unclassified-output")
                .value_name("FILE")
                .help("Filename/format to have unclassified sequences"),
        )
        .arg(
            Arg::new("kraken-output")
                .short('O')
                .long("kraken-output")
                .value_name("FILE")
                .help("Output file for normal Kraken output"),
        )
        .arg(
            Arg::new("report-kmer-data")
                .short('K')
                .long("report-kmer-data")
                .help("In comb. w/ -R, provide minimizer information in report"),
        )
        .get_matches();

    Options {
        index_filename: matches
            .get_one::<String>("index")
            .unwrap_or(&"".to_string())
            .to_string(),
        taxonomy_filename: matches
            .get_one::<String>("taxonomy")
            .unwrap_or(&"".to_string())
            .to_string(),
        options_filename: matches
            .get_one::<String>("options")
            .unwrap_or(&"".to_string())
            .to_string(),
        quick_mode: matches.contains_id("quick"),
        use_memory_mapping: matches.contains_id("memory-mapping"),
        confidence_threshold: matches
            .get_one::<String>("confidence")
            .unwrap_or(&"0".to_string())
            .parse()
            .unwrap(),
        num_threads: matches
            .get_one::<String>("threads")
            .unwrap_or(&"1".to_string())
            .parse()
            .unwrap(),
        minimum_quality_score: matches
            .get_one::<String>("min-quality")
            .unwrap_or(&"0".to_string())
            .parse()
            .unwrap(),
        paired_end_processing: matches.contains_id("paired"),
        single_file_pairs: matches.contains_id("single-file-pairs"),
        report_filename: matches
            .get_one::<String>("report")
            .unwrap_or(&"".to_string())
            .to_string(),
        mpa_style_report: matches.contains_id("mpa-style"),
        report_zero_counts: matches.contains_id("report-zero-counts"),
        print_scientific_name: matches.contains_id("scientific-name"),
        minimum_hit_groups: matches
            .get_one::<String>("hit-groups")
            .unwrap_or(&"0".to_string())
            .parse()
            .unwrap(),
        classified_output_filename: matches
            .get_one::<String>("classified-output")
            .unwrap_or(&"".to_string())
            .to_string(),
        unclassified_output_filename: matches
            .get_one::<String>("unclassified-output")
            .unwrap_or(&"".to_string())
            .to_string(),
        kraken_output_filename: matches
            .get_one::<String>("kraken-output")
            .unwrap_or(&"".to_string())
            .to_string(),
        report_kmer_data: matches.contains_id("report-kmer-data"),
        use_translated_search: false, // This will be set later based on the index options
        match_input_order: false,     // This option is not present in the command line arguments
    }
}

fn load_index_options(filename: &str) -> io::Result<IndexOptions> {
    // Implementation here
    unimplemented!();
    Ok(IndexOptions::default())
}

fn report_stats(elapsed: std::time::Duration, stats: &ClassificationStats) {
    let seconds = elapsed.as_secs_f64();
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

fn write_output(out_data: &OutputData, outputs: &Arc<Mutex<OutputStreamData>>) -> io::Result<()> {
    let outputs = outputs.lock().unwrap();
    if let Some(ref kraken_output) = outputs.kraken_output {
        kraken_output
            .lock()
            .unwrap()
            .write_all(out_data.kraken_str.as_bytes())?;
    }
    if let Some(ref classified_output1) = outputs.classified_output1 {
        classified_output1
            .lock()
            .unwrap()
            .write_all(out_data.classified_out1_str.as_bytes())?;
    }
    if let Some(ref classified_output2) = outputs.classified_output2 {
        classified_output2
            .lock()
            .unwrap()
            .write_all(out_data.classified_out2_str.as_bytes())?;
    }
    if let Some(ref unclassified_output1) = outputs.unclassified_output1 {
        unclassified_output1
            .lock()
            .unwrap()
            .write_all(out_data.unclassified_out1_str.as_bytes())?;
    }
    if let Some(ref unclassified_output2) = outputs.unclassified_output2 {
        unclassified_output2
            .lock()
            .unwrap()
            .write_all(out_data.unclassified_out2_str.as_bytes())?;
    }
    Ok(())
}

fn update_global_taxon_counters(
    thread_counters: &HashMap<TaxId, ReadCounts<readcounts::HyperLogLogPlusMinus>>,
) {
    let mut global_counters = GLOBAL_TAXON_COUNTERS.lock().unwrap();
    for (taxid, count) in thread_counters {
        global_counters
            .entry(*taxid)
            .or_insert_with(ReadCounts::default)
            .merge(count);
    }
}

fn process_sequences<'a, 'b>(
    sequences1: &'a [Sequence],
    sequences2: Option<&'a [Sequence]>,
    hash: &'a CompactHashTable,
    tax: &'a Taxonomy,
    idx_opts: &'a IndexOptions,
    opts: &'a Options,
    thread_stats: &'a mut ClassificationStats,
    scanner: &'a mut MinimizerScanner<'b>,
) -> Vec<(TaxId, String, String, String, String, String)>
where
    'a: 'b,
{
    let mut results = Vec::new();

    for (i, seq1) in sequences1.iter().enumerate() {
        let seq2 = if opts.paired_end_processing {
            // sequences2.as_ref().and_then(|seqs| seqs.get(i).cloned())
            sequences2.and_then(|seqs| seqs.get(i))
        } else {
            None
        };

        let call = classify_sequence(seq1, seq2, hash, tax, idx_opts, opts, thread_stats, scanner);

        let kraken_oss = String::new();
        let c1_oss = String::new();
        let c2_oss = String::new();
        let u1_oss = String::new();
        let u2_oss = String::new();

        // Update output strings based on classification result
        // This part would be similar to the C++ implementation
        // ...

        thread_stats.total_sequences += 1;
        thread_stats.total_bases += seq1.seq.len() as u64;
        if let Some(seq2) = seq2 {
            thread_stats.total_bases += seq2.seq.len() as u64;
        }

        results.push((call, kraken_oss, c1_oss, c2_oss, u1_oss, u2_oss));
    }

    results
}
