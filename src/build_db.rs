use std::collections::{BTreeSet, HashMap};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process;

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
    pub k: i64,
    pub l: i64,
    pub capacity: usize,
    pub maximum_capacity: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub min_clear_hash_value: u64,
    pub deterministic_build: bool,
}

// Unused in the .cc except as a forward struct, so we include it:
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

/// Placeholder for a structure that is actually defined elsewhere
/// in the C++ codebase (compact_hash.h).
#[allow(dead_code)]
struct CompactHashTable;
#[allow(dead_code)]
impl CompactHashTable {
    /// Constructor signature: `CompactHashTable(size_t capacity, int, size_t)`
    pub fn new(_capacity: usize, _some_int: i32, _bits_for_taxid: usize) -> Self {
        // no-op
        CompactHashTable
    }

    /// Stub for `WriteTable(const char*)`.
    pub fn write_table(&self, _filename: &str) {
        // no-op
    }

    /// Stub for `CompareAndSet(minimizer, new_value, &old_value)`.
    pub fn compare_and_set(&self, _key: u64, _new_val: hvalue_t, _old_val: &mut hvalue_t) -> bool {
        // always returns false just so we can mimic the loop
        false
    }

    /// Stub for `FindIndex(minimizer, &idx)`.
    pub fn find_index(&self, _key: u64, _idx: &mut usize) -> bool {
        // pretend we always fail to find
        false
    }
}

/// Placeholder for taxonomy logic (taxonomy.h).
#[allow(dead_code)]
struct Taxonomy;
#[allow(dead_code)]
impl Taxonomy {
    pub fn new(_filename: &str) -> Self {
        Taxonomy
    }
    pub fn generate_external_to_internal_id_map(&self) {
        // no-op
    }
    pub fn node_count(&self) -> usize {
        // example
        1
    }
    /// Like `taxonomy.LowestCommonAncestor(a, b)`.
    pub fn lowest_common_ancestor(&self, t1: taxid_t, t2: taxid_t) -> taxid_t {
        // trivial
        if t1 == 0 {
            t2
        } else if t2 == 0 {
            t1
        } else {
            t1.min(t2)
        }
    }
    pub fn get_internal_id(&self, ext_taxid: taxid_t) -> taxid_t {
        ext_taxid
    }
}

/// Placeholder for `NCBITaxonomy` from kraken2. In build_db.cc, used in `GenerateTaxonomy`.
#[allow(dead_code)]
struct NCBITaxonomy;
#[allow(dead_code)]
impl NCBITaxonomy {
    pub fn new(_nodes_path: &str, _names_path: &str) -> Self {
        NCBITaxonomy
    }
    pub fn mark_node(&mut self, _taxid: taxid_t) {}
    pub fn convert_to_kraken_taxonomy(&self, _filename: &str) {}
}

/// Stub for `MinimizerScanner`.
#[allow(dead_code)]
struct MinimizerScanner;
#[allow(dead_code)]
impl MinimizerScanner {
    pub fn new(_k: i64, _l: i64, _spaced_seed_mask: u64, _dna_db: bool, _toggle_mask: u64) -> Self {
        MinimizerScanner
    }
    pub fn load_sequence(&mut self, _seq: &str) {}
    pub fn load_sequence_with_range(&mut self, _seq: &str, _start: usize, _end: usize) {}
    pub fn next_minimizer(&mut self) -> Option<u64> {
        None
    }
    pub fn is_ambiguous(&self) -> bool {
        false
    }
}

/// Stub for `BatchSequenceReader`.
#[allow(dead_code)]
struct BatchSequenceReader;
#[allow(dead_code)]
impl BatchSequenceReader {
    pub fn new() -> Self {
        BatchSequenceReader
    }
    pub fn load_block(&mut self, _block_size: usize) -> bool {
        false
    }
    pub fn next_sequence(&mut self) -> Option<Sequence> {
        None
    }
}

#[allow(dead_code)]
struct Sequence {
    pub header: String,
    pub seq: String,
}

/// A placeholder for the hash function used in the .cc code (`MurmurHash3`).
#[allow(dead_code)]
fn murmur_hash3(value: u64) -> u64 {
    // dummy hash
    value ^ 0x5bd1e995
}

/// A placeholder for the ExpandSpacedSeedMask used in ParseCommandLine logic.
#[allow(dead_code)]
fn expand_spaced_seed_mask(_mask: u64, _bits_per_char: i32) {
    // no-op
}

/// Usually from “kraken2_data.h”: BITS_PER_CHAR_PRO / BITS_PER_CHAR_DNA
#[allow(dead_code)]
const BITS_PER_CHAR_PRO: i32 = 5;
#[allow(dead_code)]
const BITS_PER_CHAR_DNA: i32 = 2;

/// A placeholder for “IndexOptions” used in main() when writing out to a file.
#[allow(dead_code)]
#[repr(C)]
struct IndexOptions {
    k: i64,
    l: i64,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    dna_db: bool,
    minimum_acceptable_hash_value: u64,
    revcom_version: i32,
    db_version: i32,
    db_type: i32,
}

// Simulate these from “utilities.h” in Kraken2:
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
            let val = args[i].parse::<i64>().unwrap_or(-1);
            if val < 1 {
                errx(EX_USAGE, "k must be positive integer");
            }
            opts.k = val;
        } else if arg == "-l" {
            i += 1;
            if i >= args.len() {
                errx(EX_USAGE, "missing value for -l");
            }
            let val = args[i].parse::<i64>().unwrap_or(-1);
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
        expand_spaced_seed_mask(
            opts.spaced_seed_mask,
            if opts.input_is_protein {
                BITS_PER_CHAR_PRO
            } else {
                BITS_PER_CHAR_DNA
            },
        );
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

/// ExtractNCBISequenceIDs
#[allow(dead_code)]
fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    // from build_db.cc ~line 281
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    for &b in header.as_bytes() {
        if b == 0x01 {
            // 0x01 starts new ID
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = true;
        } else if in_id && b.is_ascii_whitespace() {
            // whitespace ends ID
            if !current_str.is_empty() {
                list.push(current_str.clone());
                current_str.clear();
            }
            in_id = false;
        } else if in_id {
            current_str.push(b as char);
        }
    }
    if !current_str.is_empty() {
        list.push(current_str);
    }
    list
}

/// SetMinimizerLCA(CompactHashTable &hash, uint64_t minimizer, taxid_t taxid, const Taxonomy &tax)
#[allow(dead_code)]
fn set_minimizer_lca(hash: &CompactHashTable, minimizer: u64, taxid: taxid_t, tax: &Taxonomy) {
    // build_db.cc line ~239
    let mut old_value: hvalue_t = 0;
    let mut new_value: hvalue_t = taxid as hvalue_t;
    // CompareAndSet loop
    while !hash.compare_and_set(minimizer, new_value, &mut old_value) {
        new_value = tax.lowest_common_ancestor(old_value as taxid_t, taxid) as hvalue_t;
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
    // build_db.cc line ~261
    scanner.load_sequence(seq);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0 && murmur_hash3(minimizer) < min_clear_hash_value {
            continue;
        }
        let mut existing_taxid = 0;
        let mut new_taxid = taxid as hvalue_t;
        // compare-and-set loop
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            new_taxid = tax.lowest_common_ancestor(existing_taxid as taxid_t, taxid) as u64;
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
        let mut block_finish = j + opts.block_size + (opts.k as usize) - 1;
        if block_finish > seq_len {
            block_finish = seq_len;
        }

        // In C++, we do subblock sets => minimizer_sets[set_ct].
        let mut minimizer_sets: Vec<BTreeSet<u64>> = (0..set_ct).map(|_| BTreeSet::new()).collect();

        // #pragma omp parallel for schedule(dynamic)
        // for (size_t i = block_start; i < block_finish; i += opts.subblock_size)
        let mut i_block = block_start;
        while i_block < block_finish {
            let subblock_finish = std::cmp::min(
                i_block + opts.subblock_size + (opts.k as usize) - 1,
                block_finish,
            );
            // local MinimizerScanner
            let mut scanner = MinimizerScanner::new(
                opts.k,
                opts.l,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
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

/// ProcessSequencesFast(...) -> ~line 148
#[allow(dead_code)]
fn process_sequences_fast(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, taxid_t>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    // #pragma omp parallel
    // We simply do a single-threaded placeholder here.
    // In real Rust code, we might use Rayon or threads.
    {
        let mut scanner = MinimizerScanner::new(
            opts.k,
            opts.l,
            opts.spaced_seed_mask,
            !opts.input_is_protein,
            opts.toggle_mask,
        );
        let mut reader = BatchSequenceReader::new();

        // while (true) ...
        loop {
            // replicate the critical LoadBlock
            let ok = reader.load_block(opts.block_size);
            if !ok {
                break;
            }
            while let Some(sequence) = reader.next_sequence() {
                let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
                let mut taxid = 0;
                for seqid in all_sequence_ids {
                    if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                        if ext_taxid != 0 {
                            taxid = taxonomy
                                .lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                        }
                    }
                }
                if taxid != 0 {
                    // Add terminator
                    let mut local_seq = sequence.seq;
                    if opts.input_is_protein && !local_seq.ends_with('*') {
                        local_seq.push('*');
                    }
                    process_sequence_fast(
                        &local_seq,
                        taxid,
                        kraken_index,
                        taxonomy,
                        &mut scanner,
                        opts.min_clear_hash_value,
                    );
                    processed_seq_ct += 1;
                    processed_ch_ct += local_seq.len();
                }
            }
        }
    }
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
}

/// ProcessSequences(...) -> ~line 185
#[allow(dead_code)]
fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, taxid_t>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    let mut reader = BatchSequenceReader::new();

    while reader.load_block(opts.block_size) {
        while let Some(mut sequence) = reader.next_sequence() {
            let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
            let mut taxid = 0;
            for seqid in all_sequence_ids {
                if let Some(&ext_taxid) = id_to_taxon_map.get(&seqid) {
                    if ext_taxid != 0 {
                        taxid = taxonomy
                            .lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                    }
                }
            }
            if taxid != 0 {
                if opts.input_is_protein && !sequence.seq.ends_with('*') {
                    sequence.seq.push('*');
                }
                process_sequence(opts, &sequence.seq, taxid, kraken_index, taxonomy);
                processed_seq_ct += 1;
                processed_ch_ct += sequence.seq.len();
            }
        }
    }
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
}

/// ReadIDToTaxonMap(...) -> ~line 429
#[allow(dead_code)]
fn read_id_to_taxon_map(id_map: &mut HashMap<String, taxid_t>, filename: &str) {
    let file = match File::open(filename) {
        Ok(f) => f,
        Err(_) => errx(EX_NOINPUT, &format!("unable to read from '{}'", filename)),
    };
    let reader = BufReader::new(file);

    for line_res in reader.lines() {
        if let Ok(line) = line_res {
            let mut parts = line.split_whitespace();
            if let Some(seq_id) = parts.next() {
                if let Some(tax_str) = parts.next() {
                    if let Ok(taxid) = tax_str.parse::<u32>() {
                        if taxid != 0 {
                            id_map.insert(seq_id.to_string(), taxid);
                        }
                    }
                }
            }
        }
    }
}

/// GenerateTaxonomy(...) -> ~line 449
#[allow(dead_code)]
fn generate_taxonomy(opts: &mut Options, id_map: &HashMap<String, taxid_t>) {
    let nodes_dmp = format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory);
    let names_dmp = format!("{}/names.dmp", opts.ncbi_taxonomy_directory);
    let mut ncbi_taxonomy = NCBITaxonomy::new(&nodes_dmp, &names_dmp);

    for (_, &taxid) in id_map.iter() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid);
        }
    }
    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename);
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
    let taxonomy = Taxonomy::new(&opts.taxonomy_filename);
    taxonomy.generate_external_to_internal_id_map();

    // bits_needed_for_value logic
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
        (32 - bits_for_taxid) as i32,
        bits_for_taxid,
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
    kraken_index.write_table(&opts.hashtable_filename);

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
    // We do a raw binary write of the struct (as in C++):
    // but that is not typically idiomatic in Rust. This is for fidelity only.
    // Warning: The layout might differ across compilers, so this is purely illustrative.
    let ptr = &index_opts as *const IndexOptions as *const u8;
    let size = std::mem::size_of::<IndexOptions>();
    let slice = unsafe { std::slice::from_raw_parts(ptr, size) };
    if writer.write_all(slice).is_err() {
        errx(EX_OSERR, "Unable to write options struct");
    }

    eprintln!("complete.");
}
