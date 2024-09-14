use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
use std::thread;

// These will likely be overridden by wrapper script,
// just keeping sane defaults in case they aren't
const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;

#[derive(Default)]
struct Options {
    id_to_taxon_map_filename: String,
    ncbi_taxonomy_directory: String,
    hashtable_filename: String,
    options_filename: String,
    taxonomy_filename: String,
    block_size: usize,
    subblock_size: usize,
    requested_bits_for_taxid: usize,
    num_threads: usize,
    input_is_protein: bool,
    k: isize,
    l: isize,
    capacity: usize,
    maximum_capacity: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    min_clear_hash_value: u64,
    deterministic_build: bool,
}

fn parse_command_line(args: Vec<String>, opts: &mut Options) {
    // parsing command line arguments here...
    // similar to the `getopt` logic in the C++ code
    // This will be implemented with a crate like `clap` or manually with `env::args()`
}

fn usage(exit_code: i32) {
    println!("Usage: build_db <options>");
    // Print other usage details...
    std::process::exit(exit_code);
}

fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = String::new();
    let mut in_id = true;

    for &ch in header.as_bytes().iter().skip(1) {
        match ch {
            0x01 => {
                if !current_str.is_empty() {
                    list.push(current_str.clone());
                }
                current_str.clear();
                in_id = true;
            }
            b' ' if in_id => {
                if !current_str.is_empty() {
                    list.push(current_str.clone());
                }
                current_str.clear();
                in_id = false;
            }
            _ if in_id => current_str.push(ch as char),
            _ => (),
        }
    }
    if !current_str.is_empty() {
        list.push(current_str);
    }

    list
}

fn set_minimizer_lca(hash: &mut CompactHashTable, minimizer: u64, taxid: taxid_t, tax: &Taxonomy) {
    let mut old_value = 0;
    let mut new_value = taxid;
    while !hash.compare_and_set(minimizer, new_value, &mut old_value) {
        new_value = tax.lowest_common_ancestor(old_value, taxid);
    }
}

fn process_sequence_fast(
    opts: &Options,
    seq: &str,
    taxid: taxid_t,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
    scanner: &mut MinimizerScanner,
) {
    scanner.load_sequence(seq);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if opts.min_clear_hash_value != 0 && murmur_hash3(minimizer) < opts.min_clear_hash_value {
            continue;
        }
        let mut existing_taxid = 0;
        let mut new_taxid = taxid;
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            new_taxid = tax.lowest_common_ancestor(new_taxid, existing_taxid);
        }
    }
}

fn process_sequences(
    opts: &Options,
    id_to_taxon_map: &HashMap<String, taxid_t>,
    hash: &mut CompactHashTable,
    tax: &Taxonomy,
) {
    let set_ct = 256;
    let locks: Vec<Mutex<()>> = (0..set_ct).map(|_| Mutex::new(())).collect();

    for j in (0..opts.capacity).step_by(opts.block_size) {
        let block_start = j;
        let block_finish = (j + opts.block_size + opts.k as usize - 1).min(opts.capacity);

        let mut minimizer_sets: Vec<HashSet<u64>> = vec![HashSet::new(); set_ct];

        (block_start..block_finish)
            .step_by(opts.subblock_size)
            .for_each(|i| {
                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    opts.spaced_seed_mask,
                    !opts.input_is_protein,
                    opts.toggle_mask,
                );
                let subblock_finish =
                    (i + opts.subblock_size + opts.k as usize - 1).min(block_finish);
                scanner.load_sequence(seq, i, subblock_finish);
                while let Some(minimizer) = scanner.next_minimizer() {
                    if scanner.is_ambiguous() {
                        continue;
                    }
                    let hc = murmur_hash3(minimizer);
                    if opts.min_clear_hash_value != 0 && hc < opts.min_clear_hash_value {
                        continue;
                    }
                    let zone = hc as usize % set_ct;
                    let _lock = locks[zone].lock().unwrap();
                    minimizer_sets[zone].insert(minimizer);
                }
            });

        let mut minimizer_list: Vec<u64> = minimizer_sets
            .into_iter()
            .flat_map(|set| set.into_iter())
            .collect();
        minimizer_list.sort_unstable();

        for &minimizer in &minimizer_list {
            set_minimizer_lca(hash, minimizer, taxid, tax);
        }
    }
}

fn read_id_to_taxon_map(id_map: &mut HashMap<String, taxid_t>, filename: &str) -> io::Result<()> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid)) = (parts.next(), parts.next()) {
            if let Ok(taxid) = taxid.parse() {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
    Ok(())
}

fn generate_taxonomy(opts: &Options, id_map: &HashMap<String, taxid_t>) {
    let mut ncbi_taxonomy = NcbiTaxonomy::new(&opts.ncbi_taxonomy_directory);
    for &taxid in id_map.values() {
        if taxid != 0 {
            ncbi_taxonomy.mark_node(taxid);
        }
    }
    ncbi_taxonomy.convert_to_kraken_taxonomy(&opts.taxonomy_filename);
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::default();
    parse_command_line(args, &mut opts);

    let id_to_taxon_map = Arc::new(Mutex::new(HashMap::new()));
    read_id_to_taxon_map(
        &mut id_to_taxon_map.lock().unwrap(),
        &opts.id_to_taxon_map_filename,
    )?;
    generate_taxonomy(&opts, &id_to_taxon_map.lock().unwrap());

    eprintln!("Taxonomy parsed and converted.");

    let taxonomy = Taxonomy::new(&opts.taxonomy_filename);
    taxonomy.generate_external_to_internal_id_map();
    let mut bits_needed_for_value = 1;
    while (1 << bits_needed_for_value) < taxonomy.node_count() as usize {
        bits_needed_for_value += 1;
    }
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        panic!("more bits required for storing taxid");
    }

    let bits_for_taxid = bits_needed_for_value.max(opts.requested_bits_for_taxid);

    let actual_capacity = if opts.maximum_capacity != 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            panic!("maximum capacity larger than requested capacity");
        }
        opts.min_clear_hash_value = ((1.0 - frac) * u64::MAX as f64) as u64;
        opts.maximum_capacity
    } else {
        opts.capacity
    };

    let mut kraken_index =
        CompactHashTable::new(actual_capacity, 32 - bits_for_taxid, bits_for_taxid);
    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    if opts.deterministic_build {
        process_sequences(
            &opts,
            &id_to_taxon_map.lock().unwrap(),
            &mut kraken_index,
            &taxonomy,
        );
    } else {
        // Assuming that the fast sequence processing function also needs to lock id_to_taxon_map.
        let id_to_taxon_map = id_to_taxon_map.lock().unwrap();
        let hash = Arc::new(Mutex::new(kraken_index));
        let tax = Arc::new(taxonomy);

        // Parallel processing using threads
        let mut handles = vec![];
        for _ in 0..opts.num_threads {
            let hash = Arc::clone(&hash);
            let tax = Arc::clone(&tax);
            let id_to_taxon_map = Arc::clone(&id_to_taxon_map);

            let handle = thread::spawn(move || {
                let mut scanner = MinimizerScanner::new(
                    opts.k,
                    opts.l,
                    opts.spaced_seed_mask,
                    !opts.input_is_protein,
                    opts.toggle_mask,
                );
                process_sequence_fast(
                    &opts,
                    &seq,
                    *taxid,
                    &mut hash.lock().unwrap(),
                    &tax,
                    &mut scanner,
                );
            });

            handles.push(handle);
        }

        for handle in handles {
            handle.join().unwrap();
        }
    }

    eprintln!("Writing data to disk...");

    kraken_index.write_table(&opts.hashtable_filename)?;

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
    let mut opts_fs = BufWriter::new(File::create(&opts.options_filename)?);
    opts_fs.write_all(bytemuck::bytes_of(&index_opts))?;

    eprintln!("Complete.");

    Ok(())
}
