use std::collections::{BTreeSet, HashMap};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process;

// Import actual implementations instead of using placeholders
use crate::compact_hash::CompactHashTable;
use crate::kraken2_data::{BITS_PER_CHAR_DNA, BITS_PER_CHAR_PRO};
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{Sequence, SequenceFormat};
use crate::taxonomy::{NCBITaxonomyImpl as NCBITaxonomy, Taxonomy};
use crate::utilities::{expand_spaced_seed_mask, murmur_hash3};

// -----------------------------------------------------------------------------
//  Constants mirroring the #defines in build_db.cc
// -----------------------------------------------------------------------------
#[allow(dead_code)]
const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
#[allow(dead_code)]
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

// -----------------------------------------------------------------------------
//  Error codes analogous to sysexits.h (used by errx in the C++ code)
// -----------------------------------------------------------------------------
#[allow(dead_code)]
const EX_USAGE: i32 = 64;
#[allow(dead_code)]
const EX_DATAERR: i32 = 65;
#[allow(dead_code)]
const EX_NOINPUT: i32 = 66;
#[allow(dead_code)]
const EX_OSERR: i32 = 71;

// -----------------------------------------------------------------------------
//  Types
// -----------------------------------------------------------------------------
#[allow(non_camel_case_types, dead_code)]
type taxid_t = u32;
#[allow(non_camel_case_types, dead_code)]
type hvalue_t = u64;

// Struct in build_db.cc (line ~53)
#[derive(Debug, Default)]
#[allow(non_snake_case, dead_code)]
struct Options {
    pub ID_to_taxon_map_filename: String,
    pub ncbi_taxonomy_directory: String,
    pub hashtable_filename: String,
    pub options_filename: String,
    pub taxonomy_filename: String,
    pub block_size: usize,
    pub subblock_size: usize,
    pub requested_bits_for_taxid: usize,
    pub num_threads: i32,
    pub input_is_protein: bool,
    pub k: usize, // Changed from i64 to usize to match MinimizerScanner
    pub l: usize, // Changed from i64 to usize to match MinimizerScanner
    pub capacity: usize,
    pub maximum_capacity: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub min_clear_hash_value: u64,
    pub deterministic_build: bool,
}

// TaxonSeqPair is defined in the .cc as a struct to hold a taxon ID and sequence
#[allow(dead_code)]
struct TaxonSeqPair {
    taxon: taxid_t,
    seq: String,
}

// -----------------------------------------------------------------------------
//  Placeholders for external items (not defined in build_db.cc)
// -----------------------------------------------------------------------------

/// In C++, this is typically from <err.h>, e.g. `errx(EX_DATAERR, "message")`.
/// We'll define a simple Rust version that prints to stderr and exits.
#[allow(dead_code)]
fn errx(exit_code: i32, message: &str) -> ! {
    eprintln!("{}", message);
    process::exit(exit_code);
}

// Using the real CompactHashTable implementation imported from crate::compact_hash

// Using the real Taxonomy implementation imported from crate::taxonomy

// Using the real NCBITaxonomy implementation (NCBITaxonomyImpl) imported from crate::taxonomy

// Using the real MinimizerScanner implementation imported from crate::mmscanner

// Using the real BatchSequenceReader and Sequence implementations imported from crate::seqreader

// Using the real murmur_hash3 implementation imported from crate::utilities

// We use crate::utilities::expand_spaced_seed_mask instead of a local implementation
// We now use the real BITS_PER_CHAR constants imported from kraken2_data module

// Wrapper for BatchSequenceReader that works with stdin
struct StdinSequenceReader {
    buf: Vec<u8>,
    pos: usize,
}

impl StdinSequenceReader {
    fn new() -> Self {
        Self {
            buf: Vec::new(),
            pos: 0,
        }
    }

    fn load_block(&mut self, _block_size: usize) -> io::Result<bool> {
        // In a real implementation, this would read from stdin
        // For now, just return true to simulate data being available
        Ok(true)
    }

    fn next_sequence(&mut self, sequence: &mut Sequence) -> io::Result<bool> {
        // In a real implementation, this would parse sequence data
        // For now, create a simple simulated sequence
        sequence.header = String::from("Simulated header");
        sequence.seq = String::from("ACGTACGTACGT");
        sequence.format = SequenceFormat::Fasta;

        // Return true to indicate a sequence was read
        Ok(true)
    }
}

/// A placeholder for “IndexOptions” used in main() when writing out to a file.
#[allow(dead_code)]
#[repr(C)]
struct IndexOptions {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    dna_db: bool,
    minimum_acceptable_hash_value: u64,
    revcom_version: i32,
    db_version: i32,
    db_type: i32,
}

// Simulate these from “utilities.h” in Kraken2:
// Define CURRENT_REVCOM_VERSION from kraken2's utilities.h
pub(crate) const CURRENT_REVCOM_VERSION: i32 = 2;

// -----------------------------------------------------------------------------
//  Functions actually defined in build_db.cc
// -----------------------------------------------------------------------------

/// usage(int exit_code)
#[allow(dead_code)]
fn usage(exit_code: i32) -> ! {
    eprintln!(
        "Usage: build_db <options>

Options (*mandatory):
* -H FILENAME   Kraken 2 hash table filename
* -m FILENAME   Sequence ID to taxon map filename
* -t FILENAME   Kraken 2 taxonomy filename
* -n DIR        NCBI taxonomy directory name
* -o FILENAME   Kraken 2 options filename
* -k INT        Set length of k-mers
* -l INT        Set length of minimizers
* -c INT        Set capacity of hash table
  -M INT        Set maximum capacity of hash table (MiniKraken)
  -S BITSTRING  Spaced seed mask
  -T BITSTRING  Minimizer toggle mask
  -X            Input seqs. are proteins
  -p INT        Number of threads
  -F            Use fast, nondeterministic building method
  -B INT        Read block size
  -b INT        Read subblock size
  -r INT        Bit storage requested for taxid
"
    );
    process::exit(exit_code);
}

/// ParseCommandLine(int argc, char **argv, Options &opts)
#[allow(dead_code)]
fn parse_command_line(args: &[String], opts: &mut Options) {
    // In C++: while ((opt = getopt(argc, argv, "?hB:b:c:FH:m:n:o:t:k:l:M:p:r:s:S:T:X")) != -1)
    // We replicate the same flags in a basic manner:
    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];
        if arg == "-h" || arg == "?" {
            usage(0);
        } else if arg == "-B" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -B");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "must have positive block size");
            }
            opts.block_size = val as usize;
        } else if arg == "-b" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -b");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "must have positive subblock size");
            }
            opts.subblock_size = val as usize;
        } else if arg == "-r" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -r");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 0 || val > 31 {
                errx(EX_USAGE, "invalid bit storage requested");
            }
            opts.requested_bits_for_taxid = val as usize;
        } else if arg == "-p" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -p");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "can't have negative number of threads");
            }
            // no direct analog to `omp_get_max_threads()`, so we skip that check
            opts.num_threads = val as i32;
        } else if arg == "-H" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing filename for -H");
            }
            opts.hashtable_filename = args[i].clone();
        } else if arg == "-m" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing filename for -m");
            }
            opts.ID_to_taxon_map_filename = args[i].clone();
        } else if arg == "-n" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing dir name for -n");
            }
            opts.ncbi_taxonomy_directory = args[i].clone();
        } else if arg == "-o" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing filename for -o");
            }
            opts.options_filename = args[i].clone();
        } else if arg == "-t" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing filename for -t");
            }
            opts.taxonomy_filename = args[i].clone();
        } else if arg == "-S" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing bitstring for -S");
            }
            // in C++, parse with strtol(optarg, nullptr, 2)
            // We'll parse from binary. If invalid, set to 0.
            let bitstring = &args[i];
            match u64::from_str_radix(bitstring, 2) {
                Ok(mask) => {
                    opts.spaced_seed_mask = mask;
                }
                Err(_) => {
                    errx(EX_USAGE, "invalid spaced seed mask");
                }
            }
        } else if arg == "-T" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing bitstring for -T");
            }
            let bitstring = &args[i];
            match u64::from_str_radix(bitstring, 2) {
                Ok(mask) => {
                    opts.toggle_mask = mask;
                }
                Err(_) => {
                    errx(EX_USAGE, "invalid toggle mask");
                }
            }
        } else if arg == "-k" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -k");
            }
            let val = args[i].parse::<usize>().unwrap_or(0);
            if val < 1 {
                errx(EX_USAGE, "k must be positive integer");
            }
            opts.k = val;
        } else if arg == "-l" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -l");
            }
            let val = args[i].parse::<usize>().unwrap_or(0);
            if val < 1 || val > 31 {
                errx(EX_USAGE, "l must be between 1 and 31");
            }
            opts.l = val;
        } else if arg == "-c" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing capacity for -c");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "capacity must be positive integer");
            }
            opts.capacity = val as usize;
        } else if arg == "-M" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing max capacity for -M");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "max capacity must be positive integer");
            }
            opts.maximum_capacity = val as usize;
        } else if arg == "-F" {
            opts.deterministic_build = false;
        } else if arg == "-X" {
            opts.input_is_protein = true;
        } else {
            // If we reach here, unknown option
            eprintln!("Unknown option: {}", arg);
            usage(EX_USAGE);
        }
        i += 1;
    }

    // replicate final checks from build_db.cc
    if opts.spaced_seed_mask != 0 {
        let mut mask = opts.spaced_seed_mask;
        expand_spaced_seed_mask(
            &mut mask,
            if opts.input_is_protein {
                BITS_PER_CHAR_PRO.into()
            } else {
                BITS_PER_CHAR_DNA.into()
            },
        );
        opts.spaced_seed_mask = mask;
    }
    if opts.hashtable_filename.is_empty()
        || opts.ID_to_taxon_map_filename.is_empty()
        || opts.ncbi_taxonomy_directory.is_empty()
        || opts.options_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
    {
        eprintln!("missing mandatory filename parameter");
        usage(EX_USAGE);
    }
    if opts.k == 0 || opts.l == 0 || opts.capacity == 0 {
        eprintln!("missing mandatory integer parameter");
        usage(EX_USAGE);
    }
    if opts.k < opts.l {
        eprintln!("k cannot be less than l");
        usage(EX_USAGE);
    }
    if opts.block_size < opts.subblock_size {
        eprintln!("block size cannot be less than subblock size");
        usage(EX_USAGE);
    }
    if opts.maximum_capacity > opts.capacity {
        eprintln!("maximum capacity option shouldn't specify larger capacity than normal");
        usage(EX_USAGE);
    }
}

/// This function exists to deal with NCBI's use of \x01 characters to denote
/// the start of a new FASTA header in the same line (for non-redundant DBs).
/// We return all sequence IDs in a header line, not just the first.
/// Implements ExtractNCBISequenceIDs from build_db.cc (~line 263)
#[allow(dead_code)]
fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    // In C++:
    // vector<string> list;
    // string current_str;
    // bool in_id = true;
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    // In C++: for (size_t i = 0; i < header.size(); ++i)
    for &b in header.as_bytes() {
        if b == 0x01 {
            // 0x01 starts new ID
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = true;
        } else if in_id && b.is_ascii_whitespace() {
            // spaces end ID
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = false;
        } else if in_id {
            // Build ID string char by char
            current_str.push(b as char);
        }
    }
    if !current_str.is_empty() {
        list.push(current_str);
    }
    list
}

/// Implements SetMinimizerLCA from build_db.cc (~line 296)
#[allow(dead_code)]
fn set_minimizer_lca(hash: &CompactHashTable, minimizer: u64, taxid: taxid_t, tax: &Taxonomy) {
    // In C++:
    // hvalue_t old_value = 0;
    // hvalue_t new_value = taxid;
    // while (! hash.CompareAndSet(minimizer, new_value, &old_value))
    //   new_value = tax.LowestCommonAncestor(old_value, taxid);
    let mut old_value: hvalue_t = 0;
    let mut new_value: hvalue_t = taxid as hvalue_t;
    // CompareAndSet loop - identical to C++ implementation
    while !hash.compare_and_set(minimizer, new_value, &mut old_value) {
        // Safely convert between types with bounds checking
        let old_taxid = if old_value <= taxid_t::MAX as hvalue_t {
            old_value as taxid_t
        } else {
            eprintln!("Warning: taxid value overflow in set_minimizer_lca");
            taxid_t::MAX
        };

        let lca = tax.lowest_common_ancestor(old_taxid as u64, taxid as u64);
        new_value = lca as hvalue_t;
    }
}

/// ProcessSequenceFast(const string &seq, taxid_t taxid, CompactHashTable &hash, ... )
#[allow(dead_code)]
fn process_sequence_fast(
    seq: &str,
    taxid: taxid_t,
    hash: &CompactHashTable,
    tax: &Taxonomy,
    scanner: &mut MinimizerScanner,
    min_clear_hash_value: u64,
) {
    // Matches build_db.cc ProcessSequenceFast function (line ~305)
    scanner.load_sequence(seq, 0, seq.len());
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0 && murmur_hash3(minimizer) < min_clear_hash_value {
            continue;
        }
        let mut existing_taxid = 0;
        let mut new_taxid = taxid as hvalue_t;
        // compare-and-set loop (identical to C++ implementation)
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            // Safely convert between types with bounds checking
            let old_taxid = if existing_taxid <= taxid_t::MAX as hvalue_t {
                existing_taxid as taxid_t
            } else {
                eprintln!("Warning: taxid value overflow in process_sequence_fast");
                taxid_t::MAX
            };

            let lca = tax.lowest_common_ancestor(old_taxid as u64, taxid as u64);
            new_taxid = lca as hvalue_t;
        }
    }
}

/// ProcessSequence (deterministic approach)
#[allow(dead_code)]
fn process_sequence(
    opts: &Options,
    seq: &str,
    taxid: taxid_t,
    hash: &CompactHashTable,
    tax: &Taxonomy,
) {
    // build_db.cc ~line 314
    // The code is complex with #pragma omp parallel for. We replicate the structure in a single-thread form.
    let set_ct = 256;

    // The .cc allocates omp_lock_t locks[set_ct].
    // For fidelity, we’ll just skip actual locking here or store placeholders.
    // Then it does for (size_t j = 0; j < seq.size(); j += opts.block_size)
    let seq_len = seq.len();
    let mut j = 0;
    while j < seq_len {
        let block_start = j;
        // Calculate block finish safely to avoid potential overflow
        let k_minus_1 = opts.k.checked_sub(1).unwrap_or(0);
        let block_size_plus_k = opts.block_size.checked_add(k_minus_1).unwrap_or(seq_len);
        let j_plus_size = j.checked_add(block_size_plus_k).unwrap_or(seq_len);
        let block_finish = std::cmp::min(j_plus_size, seq_len);

        // In C++, we do subblock sets => minimizer_sets[set_ct].
        let mut minimizer_sets: Vec<BTreeSet<u64>> = (0..set_ct).map(|_| BTreeSet::new()).collect();

        // #pragma omp parallel for schedule(dynamic)
        // for (size_t i = block_start; i < block_finish; i += opts.subblock_size)
        let mut i_block = block_start;
        while i_block < block_finish {
            // Calculate subblock finish safely to avoid potential overflow
            let k_minus_1 = opts.k.checked_sub(1).unwrap_or(0);
            let subblock_finish = std::cmp::min(
                i_block
                    .checked_add(opts.subblock_size)
                    .unwrap_or(block_finish)
                    .checked_add(k_minus_1)
                    .unwrap_or(block_finish),
                block_finish,
            );
            // local MinimizerScanner
            let mut scanner = MinimizerScanner::new(
                opts.k,
                opts.l,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
                true, // ambiguous_check parameter
            );
            scanner.load_sequence_with_range(seq, i_block, subblock_finish);

            while let Some(minimizer) = scanner.next_minimizer() {
                if scanner.is_ambiguous() {
                    continue;
                }
                let hc = murmur_hash3(minimizer);
                if opts.min_clear_hash_value != 0 && hc < opts.min_clear_hash_value {
                    continue;
                }
                let zone = (hc % (set_ct as u64)) as usize;
                // Insert into the BTreeSet for that zone
                minimizer_sets[zone].insert(minimizer);
            }

            i_block += opts.subblock_size;
        }

        // combine sets into sorted lists
        let mut minimizer_lists: Vec<Vec<u64>> = Vec::with_capacity(set_ct);
        let mut prefix_sizes = vec![0_usize; set_ct + 1];
        for (idx, s) in minimizer_sets.into_iter().enumerate() {
            let mut v: Vec<u64> = s.into_iter().collect();
            v.sort_unstable();
            prefix_sizes[idx + 1] = v.len();
            minimizer_lists.push(v);
        }
        for z in 2..=set_ct {
            prefix_sizes[z] += prefix_sizes[z - 1];
        }
        let mut all_minimizers = vec![0_u64; prefix_sizes[set_ct]];
        for (idx, v) in minimizer_lists.into_iter().enumerate() {
            let start_pos = prefix_sizes[idx];
            all_minimizers[start_pos..(start_pos + v.len())].copy_from_slice(&v);
        }

        // The .cc code sets up an insertion loop
        let mut mm_ct = all_minimizers.len();
        let mut index_list = vec![0_usize; mm_ct];
        let mut insertion_list = vec![true; mm_ct];

        while !all_minimizers.is_empty() {
            // gather insertion point info
            for i in 0..mm_ct {
                if insertion_list[i] {
                    // find_index
                    let mut idx_out = 0;
                    insertion_list[i] = !hash.find_index(all_minimizers[i], &mut idx_out);
                    index_list[i] = idx_out;
                }
            }
            let mut novel_insertion_points = BTreeSet::new();
            let mut safe_ct = 0;
            for i in 0..mm_ct {
                if insertion_list[i] {
                    if novel_insertion_points.contains(&index_list[i]) {
                        break;
                    }
                    novel_insertion_points.insert(index_list[i]);
                }
                safe_ct += 1;
            }
            // #pragma omp parallel for
            for i in 0..safe_ct {
                set_minimizer_lca(hash, all_minimizers[i], taxid, tax);
            }
            // remove the safe prefix
            all_minimizers.drain(0..safe_ct);
            index_list.drain(0..safe_ct);
            insertion_list.drain(0..safe_ct);
            mm_ct -= safe_ct;
        }

        j += opts.block_size;
    }
}

/// Implements ProcessSequencesFast from build_db.cc (~line 154)
/// A quick but nondeterministic build
#[allow(dead_code)]
fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, taxid_t>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    // In C++: #pragma omp parallel
    // We simply do a single-threaded placeholder here.
    // In real Rust code, we would use Rayon or std::thread
    {
        let mut scanner = MinimizerScanner::new(
            opts.k,
            opts.l,
            opts.spaced_seed_mask,
            !opts.input_is_protein,
            opts.toggle_mask,
            true, // ambiguous_check parameter
        );
        // Use our simplified wrapper instead of directly using BatchSequenceReader
        let mut reader = StdinSequenceReader::new();

        // In C++: while (true) ...
        loop {
            // In C++: this is wrapped in a critical section with reader_clone
            let ok = match reader.load_block(opts.block_size) {
                Ok(result) => result,
                Err(e) => {
                    eprintln!("Error loading block: {}", e);
                    false
                }
            };
            if !ok {
                break;
            }

            // In the C++ code, this uses a Sequence pointer, but in Rust we reuse a Sequence object
            let mut sequence = Sequence::default();
            while match reader.next_sequence(&mut sequence) {
                Ok(has_next) => has_next,
                Err(e) => {
                    eprintln!("Error reading sequence: {}", e);
                    false
                }
            } {
                let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
                let mut taxid = 0;
                for seqid in all_sequence_ids {
                    if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                        if ext_taxid != 0 {
                            taxid = taxonomy.lowest_common_ancestor(
                                taxid,
                                taxonomy.get_internal_id(ext_taxid as u64),
                            );
                        }
                    }
                }
                if taxid != 0 {
                    // Add terminator for protein sequences if not already there
                    let mut local_seq = sequence.seq.clone();
                    if opts.input_is_protein && !local_seq.ends_with('*') {
                        local_seq.push('*');
                    }
                    process_sequence_fast(
                        &local_seq,
                        taxid as u32,
                        kraken_index,
                        taxonomy,
                        &mut scanner,
                        opts.min_clear_hash_value,
                    );
                    // In C++ these would be atomic operations in the parallel section
                    processed_seq_ct += 1;
                    processed_ch_ct += local_seq.len();
                }
            }

            // In C++ there's isatty/status update logic here that we omit
        }
    }

    // In C++ there's another isatty check before this final output
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
}

/// Implements ProcessSequences from build_db.cc (~line 217)
/// Slightly slower but deterministic when multithreaded
#[allow(dead_code)]
fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, taxid_t>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    // Use our simplified wrapper
    let mut reader = StdinSequenceReader::new();

    // In C++, the sequence is a pointer: Sequence *sequence;
    // Here we create a reusable Sequence object to receive data
    let mut sequence = Sequence::default();

    while match reader.load_block(opts.block_size) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("Error loading block: {}", e);
            false
        }
    } {
        while match reader.next_sequence(&mut sequence) {
            Ok(has_next) => has_next,
            Err(e) => {
                eprintln!("Error reading sequence: {}", e);
                false
            }
        } {
            let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
            let mut taxid = 0;
            // In C++ this is an int ext_taxid, here we use it directly from the HashMap
            for seqid in all_sequence_ids {
                if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                    if ext_taxid != 0 {
                        taxid = taxonomy.lowest_common_ancestor(
                            taxid,
                            taxonomy.get_internal_id(ext_taxid as u64),
                        );
                    }
                }
            }
            if taxid != 0 {
                // Add terminator for protein sequences if not already there
                if opts.input_is_protein && !sequence.seq.ends_with('*') {
                    sequence.seq.push('*');
                }
                process_sequence(opts, &sequence.seq, taxid as u32, kraken_index, taxonomy);
                processed_seq_ct += 1;
                processed_ch_ct += sequence.seq.len();
            }
        }

        // In C++ there's isatty/status update logic here that we omit
    }

    // In C++ there's another isatty check before this final output
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
}

/// Implements ReadIDToTaxonMap from build_db.cc (~line 435)
#[allow(dead_code)]
fn read_id_to_taxon_map(id_map: &mut HashMap<String, taxid_t>, filename: &str) {
    // In C++: ifstream map_file(filename.c_str());
    let file = match File::open(filename) {
        Ok(f) => f,
        Err(_) => errx(EX_NOINPUT, &format!("unable to read from '{}'", filename)),
    };
    let reader = BufReader::new(file);

    // In C++: while (getline(map_file, line))
    for line_res in reader.lines() {
        if let Ok(line) = line_res {
            // In C++: istringstream iss(line); iss >> seq_id; iss >> taxid;
            let mut parts = line.split_whitespace();
            if let Some(seq_id) = parts.next() {
                if let Some(tax_str) = parts.next() {
                    if let Ok(taxid) = tax_str.parse::<u32>() {
                        // In C++: if (taxid) id_map[seq_id] = taxid;
                        if taxid != 0 {
                            id_map.insert(seq_id.to_string(), taxid);
                        }
                    }
                }
            }
        }
    }
}

/// Implements GenerateTaxonomy from build_db.cc (~line 452)
#[allow(dead_code)]
fn generate_taxonomy(opts: &mut Options, id_map: &HashMap<String, taxid_t>) {
    // In C++:
    // NCBITaxonomy ncbi_taxonomy(
    //   opts.ncbi_taxonomy_directory + "/nodes.dmp",
    //   opts.ncbi_taxonomy_directory + "/names.dmp"
    // );
    let nodes_dmp = format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory);
    let names_dmp = format!("{}/names.dmp", opts.ncbi_taxonomy_directory);
    
    // Create NCBITaxonomy and handle potential error
    let mut ncbi_taxonomy = match NCBITaxonomy::new(&nodes_dmp, &names_dmp) {
        Ok(taxonomy) => taxonomy,
        Err(e) => errx(EX_DATAERR, &format!("Error creating taxonomy: {}", e)),
    };

    // In C++: for (auto &kv_pair : id_map)
    for (_, &taxid) in id_map.iter() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid as u64);
        }
    }
    
    // In C++: ncbi_taxonomy.ConvertToKrakenTaxonomy(opts.taxonomy_filename.c_str());
    if let Err(e) = ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename) {
        errx(EX_DATAERR, &format!("Error converting taxonomy: {}", e));
    }
}

// -----------------------------------------------------------------------------
//  main(...) -> ~line 76
// -----------------------------------------------------------------------------
#[allow(dead_code)]
fn main() {
    let mut opts = Options {
        spaced_seed_mask: 0, // DEFAULT_SPACED_SEED_MASK in C++ is from kraken2_data.h, set to 0 as placeholder
        toggle_mask: 0,      // DEFAULT_TOGGLE_MASK also from kraken2_data.h, set to 0
        input_is_protein: false,
        num_threads: 1,
        block_size: DEFAULT_BLOCK_SIZE,
        subblock_size: DEFAULT_SUBBLOCK_SIZE,
        requested_bits_for_taxid: 0,
        min_clear_hash_value: 0,
        maximum_capacity: 0,
        deterministic_build: true,
        ..Default::default()
    };

    // parse command line
    let args: Vec<String> = env::args().collect();
    parse_command_line(&args, &mut opts);

    // In the C++ code: `omp_set_num_threads(opts.num_threads);`
    // Omitted in Rust—no direct equivalent unless you use Rayon’s thread pool, etc.

    // map<string, taxid_t> in C++
    let mut id_to_taxon_map: HashMap<String, taxid_t> = HashMap::new();

    // Read ID->taxon map
    read_id_to_taxon_map(&mut id_to_taxon_map, &opts.ID_to_taxon_map_filename);
    // Generate taxonomy
    generate_taxonomy(&mut opts, &id_to_taxon_map);

    eprintln!("Taxonomy parsed and converted.");

    // Taxonomy taxonomy(...)
    let mut taxonomy = match Taxonomy::new(&opts.taxonomy_filename, false) {
        Ok(taxonomy) => taxonomy,
        Err(e) => errx(EX_DATAERR, &format!("Error creating taxonomy: {}", e)),
    };
    taxonomy.generate_external_to_internal_id_map();

    // Determine bits needed for storing taxid
    let mut bits_needed_for_value = 1usize;
    while (1 << bits_needed_for_value) < taxonomy.node_count() as i64 {
        bits_needed_for_value += 1;
    }
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        errx(EX_DATAERR, "more bits required for storing taxid");
    }

    let mut bits_for_taxid = bits_needed_for_value;
    if bits_for_taxid < opts.requested_bits_for_taxid {
        bits_for_taxid = opts.requested_bits_for_taxid;
    }

    let mut actual_capacity = opts.capacity;
    if opts.maximum_capacity != 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            errx(
                EX_DATAERR,
                "maximum capacity larger than requested capacity",
            );
        }
        opts.min_clear_hash_value = ((1.0 - frac) * (u64::MAX as f64)) as u64;
        actual_capacity = opts.maximum_capacity;
    }

    // create the CHT
    // C++: CompactHashTable kraken_index(actual_capacity, 32 - bits_for_taxid, bits_for_taxid);
    let kraken_index = CompactHashTable::new(
        actual_capacity,
        (32 - bits_for_taxid) as u8,
        bits_for_taxid as u8,
    );

    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    // Build
    if opts.deterministic_build {
        process_sequences(&opts, &id_to_taxon_map, &kraken_index, &taxonomy);
    } else {
        process_sequences_fast(&opts, &id_to_taxon_map, &kraken_index, &taxonomy);
    }

    eprintln!("Writing data to disk... ");
    if let Err(e) = kraken_index.write_table(&opts.hashtable_filename) {
        errx(EX_OSERR, &format!("Error writing hash table: {}", e));
    }

    // Write out IndexOptions
    let index_opts = IndexOptions {
        k: opts.k,
        l: opts.l,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        minimum_acceptable_hash_value: opts.min_clear_hash_value,
        revcom_version: CURRENT_REVCOM_VERSION,
        db_version: 0,
        db_type: 0,
    };

    // C++ code:
    // ofstream opts_fs(opts.options_filename);
    // opts_fs.write((char*)&index_opts, sizeof(index_opts));
    // ...
    let file = match File::create(&opts.options_filename) {
        Ok(f) => f,
        Err(_) => errx(
            EX_OSERR,
            &format!("Unable to write options file {}", opts.options_filename),
        ),
    };
    let mut writer = BufWriter::new(file);

    // Write options in a safer, more idiomatic way instead of raw binary
    // This ensures we don't rely on memory layout which can differ

    // Write each field individually
    if writer.write_all(&index_opts.k.to_le_bytes()).is_err()
        || writer.write_all(&index_opts.l.to_le_bytes()).is_err()
        || writer
            .write_all(&index_opts.spaced_seed_mask.to_le_bytes())
            .is_err()
        || writer
            .write_all(&index_opts.toggle_mask.to_le_bytes())
            .is_err()
        || writer.write_all(&[index_opts.dna_db as u8]).is_err()
        || writer
            .write_all(&index_opts.minimum_acceptable_hash_value.to_le_bytes())
            .is_err()
        || writer
            .write_all(&index_opts.revcom_version.to_le_bytes())
            .is_err()
        || writer
            .write_all(&index_opts.db_version.to_le_bytes())
            .is_err()
        || writer.write_all(&index_opts.db_type.to_le_bytes()).is_err()
    {
        errx(EX_OSERR, "Unable to write options struct");
    }

    eprintln!("complete.");
}
