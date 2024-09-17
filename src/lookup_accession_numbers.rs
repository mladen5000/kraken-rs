/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process;

fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lookup_accession_numbers <lookup file> <accmaps>");
        process::exit(1);
    }

    // Read the lookup list file
    let lookup_file = &args[1];
    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();
    let lookup_list_file = File::open(lookup_file).unwrap_or_else(|_| {
        eprintln!("Error opening lookup file: {}", lookup_file);
        process::exit(1);
    });

    let reader = BufReader::new(lookup_list_file);
    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.splitn(2, '\t').collect();
        if fields.len() != 2 {
            continue;
        }
        let seqid = fields[0].to_string();
        let accnum = fields[1].to_string();
        target_lists
            .entry(accnum)
            .or_default()
            .push(seqid);
    }

    let initial_target_count = target_lists.len();
    let mut accessions_searched: u64 = 0;

    if atty::is(atty::Stream::Stderr) {
        eprint!("\rFound 0/{} targets...", initial_target_count);
    }

    // Process each accmap file
    for accmap_filename in &args[2..] {
        if target_lists.is_empty() {
            // Stop processing files if we've found all we need
            break;
        }

        let accmap_path = Path::new(accmap_filename);
        let accmap_file = File::open(accmap_path).unwrap_or_else(|_| {
            eprintln!("Error opening accmap file: {}", accmap_filename);
            process::exit(1);
        });

        let metadata = accmap_file.metadata().unwrap();
        let filesize = metadata.len();
        let mmap = unsafe {
            memmap2::Mmap::map(&accmap_file).unwrap_or_else(|_| {
                eprintln!("Error memory-mapping file: {}", accmap_filename);
                process::exit(1);
            })
        };

        let mut ptr = &mmap[..];

        // Skip header line
        if let Some(pos) = ptr.iter().position(|&c| c == b'\n') {
            ptr = &ptr[pos + 1..];
        }

        while !ptr.is_empty() {
            if target_lists.is_empty() {
                // Stop processing file if we've found all we need
                break;
            }

            // Find the end of the current line
            let line_end = match ptr.iter().position(|&c| c == b'\n') {
                Some(pos) => pos,
                None => {
                    eprintln!("expected EOL not found at EOF in {}", accmap_filename);
                    break;
                }
            };

            // Get the line
            let line = &ptr[..line_end];
            ptr = &ptr[line_end + 1..];

            // Split the line into fields
            let mut fields_iter = line.split(|&c| c == b'\t');
            let accnum_field = fields_iter.next();
            if accnum_field.is_none() {
                eprintln!("expected TAB not found in {}", accmap_filename);
                break;
            }
            let accnum = String::from_utf8_lossy(accnum_field.unwrap()).to_string();
            accessions_searched += 1;

            if target_lists.contains_key(&accnum) {
                // Skip the next two fields to get to the taxid
                let mut taxid_field = None;
                for _ in 0..2 {
                    taxid_field = fields_iter.next();
                    if taxid_field.is_none() {
                        eprintln!("expected TAB not found in {}", accmap_filename);
                        break;
                    }
                }
                if taxid_field.is_none() {
                    break;
                }
                let taxid = String::from_utf8_lossy(taxid_field.unwrap()).to_string();
                let seqids = target_lists.remove(&accnum).unwrap();
                for seqid in seqids {
                    println!("{}\t{}", seqid, taxid);
                }
                if atty::is(atty::Stream::Stderr) {
                    eprint!(
                        "\rFound {}/{} targets, searched through {} accession IDs...",
                        initial_target_count - target_lists.len(),
                        initial_target_count,
                        accessions_searched
                    );
                }
            }

            if accessions_searched % 10_000_000 == 0 && atty::is(atty::Stream::Stderr) {
                eprint!(
                    "\rFound {}/{} targets, searched through {} accession IDs...",
                    initial_target_count - target_lists.len(),
                    initial_target_count,
                    accessions_searched
                );
            }
        }
    }

    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }
    eprintln!(
        "Found {}/{} targets, searched through {} accession IDs, search complete.",
        initial_target_count - target_lists.len(),
        initial_target_count,
        accessions_searched
    );

    if !target_lists.is_empty() {
        eprintln!(
            "lookup_accession_numbers: {}/{} accession numbers remain unmapped, see \
              unmapped.txt in DB directory",
            target_lists.len(),
            initial_target_count
        );
        let mut ofs = File::create("unmapped.txt").unwrap();
        for accnum in target_lists.keys() {
            writeln!(ofs, "{}", accnum).unwrap();
        }
    }
}
