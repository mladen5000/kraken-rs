use crate::ncbi_taxonomy;
use crate::taxonomy::NCBITaxonomy;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;
const DEFAULT_BLOCK_SIZE: usize = (10 * 1024 * 1024);
const DEFAULT_SUBBLOCK_SIZE: usize = (1024);

struct Options {
    ID_to_taxon_map_filename: String,
    ncbi_taxonomy_directory: String,
    hashtable_filename: String,
    options_filename: String,
    taxonomy_filename: String,
    block_size: usize,
    subblock_size: usize,
    k: isize,
    l: isize,
    requested_bits_for_taxid: usize,
    num_threads: i64,
    input_is_protein: bool,
    capacity: usize,
    maximum_capacity: usize,
    spaced_seed_mask: i64,
    toggle_mask: i64,
    min_clear_hash_value: i64,
    deterministic_build: bool,
}

type TaxidT = usize;

fn read_id_to_taxon_map(id_map: &mut HashMap<String, TaxidT>, filename: &str) {
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(_) => {
            eprintln!("Unable to read from '{}'", filename);
            process::exit(1);
        }
    };

    let reader = BufReader::new(file);

    // Loop through each line
    for line_result in reader.lines() {
        let line = line_result.expect("Error reading line");
        let mut parts = line.split_whitespace();
        if let (Some(seq_id), Some(taxid_str)) = (parts.next(), parts.next()) {
            if let Ok(taxid) = taxid_str.parse::<TaxidT>() {
                id_map.insert(seq_id.to_string(), taxid);
            }
        }
    }
}

/*
fn generate_taxonomy(opts: Options, id_map: &mut HashMap<String, TaxidT>) {
    let ncbi_taxonomy = NCBITaxonomy::new(
        &format!("{}/nodes.dmp", opts.ncbi_taxonomy_directory),
        &format!("{}/names.dmp", opts.ncbi_taxonomy_directory),
    );

    for kv_pair in id_map {
        if *kv_pair.1 != 0 {
            ncbi_taxonomy.mark_node(*kv_pair.1);
        }
    }
}
*/
