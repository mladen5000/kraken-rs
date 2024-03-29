/*
* Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
*
* This file is part of the Kraken 2 taxonomic sequence classification system.
*/

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lookup_accession_numbers <lookup file> <accmaps>");
        std::process::exit(1);
    }

    //Map from accession number to list of sequence IDs
    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();
    let lookup_file = File::open(&args[1])?;
    let reader = BufReader::new(lookup_file);
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let seqid = fields[0].to_string();
        let accnum = fields[1].to_string();
        target_lists.entry(accnum).or_default().push(seqid);
    }

    let initial_target_count = target_lists.len();

    let mut accessions_searched = 0;
    // Iterate over the accession maps
    for accmap in &args[2..] {
        if target_lists.is_empty() {
            break;
        }

        // Open the accession map file and search for accession numbers and their corresponding taxids
        let accmap_file = File::open(accmap)?;
        let reader = BufReader::new(accmap_file);

        // Iterate over the lines in the accession map file
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            let accnum = fields[0].to_string();
            accessions_searched += 1;

            if let Some(seqids) = target_lists.remove(&accnum) {
                let taxid = fields[2].to_string();
                for seqid in seqids {
                    println!("{}\t{}", seqid, taxid);
                }
            }

            if accessions_searched % 10_000_000 == 0 {
                eprintln!(
                    "\rFound {}/{} targets, searched through {} accession IDs...",
                    initial_target_count - target_lists.len(),
                    initial_target_count,
                    accessions_searched
                );
            }
        }
    }

    eprintln!(
        "Found {}/{} targets, searched through {} accession IDs, search complete.",
        initial_target_count - target_lists.len(),
        initial_target_count,
        accessions_searched
    );

    if !target_lists.is_empty() {
        eprintln!(
            "lookup_accession_numbers: {}/{} accession numbers remain unmapped, see unmapped.txt in DB directory",
            target_lists.len(),
            initial_target_count
        );
        let mut file = File::create("unmapped.txt")?;
        for accnum in target_lists.keys() {
            writeln!(file, "{}", accnum)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_lookup_accession_numbers() {
        // Create test input files
        let lookup_file_content = "seq1\tacc1\nseq2\tacc2\nseq3\tacc1\n";
        let accmap_file_content = "acc1\t123\t1000\nacc2\t456\t2000\n";

        let lookup_file_path = "test_lookup.txt";
        let accmap_file_path = "test_accmap.txt";

        let mut lookup_file = File::create(lookup_file_path).unwrap();
        lookup_file
            .write_all(lookup_file_content.as_bytes())
            .unwrap();

        let mut accmap_file = File::create(accmap_file_path).unwrap();
        accmap_file
            .write_all(accmap_file_content.as_bytes())
            .unwrap();

        // Call the main function with test arguments
        let args = vec![
            "lookup_accession_numbers".to_string(),
            lookup_file_path.to_string(),
            accmap_file_path.to_string(),
        ];
        std::env::set_var("RUST_BACKTRACE", "1");
        main().unwrap();

        // Check the output
        let expected_output: &str = "seq1\t1000\nseq3\t1000\nseq2\t2000\n";
        assert_eq!(std::fs::read_to_string("unmapped.txt").unwrap(), "");

        // Clean up test files
        std::fs::remove_file(lookup_file_path).unwrap();
        std::fs::remove_file(accmap_file_path).unwrap();
        std::fs::remove_file("unmapped.txt").unwrap();
    }

    #[test]
    fn test_lookup_accession_numbers_with_unmapped() {
        // Create test input files
        let lookup_file_content = "seq1\tacc1\nseq2\tacc2\nseq3\tacc3\n";
        let accmap_file_content = "acc1\t123\t1000\n";

        let lookup_file_path = "test_lookup.txt";
        let accmap_file_path = "test_accmap.txt";

        let mut lookup_file = File::create(lookup_file_path).unwrap();
        lookup_file
            .write_all(lookup_file_content.as_bytes())
            .unwrap();

        let mut accmap_file = File::create(accmap_file_path).unwrap();
        accmap_file
            .write_all(accmap_file_content.as_bytes())
            .unwrap();

        // Call the main function with test arguments
        let args = vec![
            "lookup_accession_numbers".to_string(),
            lookup_file_path.to_string(),
            accmap_file_path.to_string(),
        ];
        std::env::set_var("RUST_BACKTRACE", "1");
        main().unwrap();

        // Check the output
        let expected_output = "seq1\t1000\n";
        let expected_unmapped = "acc2\nacc3\n";
        assert_eq!(
            std::fs::read_to_string("unmapped.txt").unwrap(),
            expected_unmapped
        );

        // Clean up test files
        std::fs::remove_file(lookup_file_path).unwrap();
        std::fs::remove_file(accmap_file_path).unwrap();
        std::fs::remove_file("unmapped.txt").unwrap();
    }
}
