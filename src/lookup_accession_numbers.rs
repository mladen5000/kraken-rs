/*
* Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
*
* This file is part of the Kraken 2 taxonomic sequence classification system.
*/

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Main functions performs the following:
/// 1. Read the lookup file and create a map from accession number to list of sequence IDs
/// 2. Iterate over the accession maps and search for accession numbers and their corresponding taxids
/// 3. Print the sequence ID and taxid to stdout
/// 4. Print the unmapped accession numbers to a file
pub fn main() -> io::Result<()> {
    // Gather args (lookupfile) and (accession maps)
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lookup_accession_numbers <lookup file> <accmaps>");
        std::process::exit(1);
    }

    // Open lookup file
    // Map from accession number to list of sequence IDs
    let lookup_file = File::open(&args[1])?;
    let reader = BufReader::new(lookup_file);

    // Get seqID and AccessionNUM from tsv data
    // Build a accession to sequence ID table
    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let seqid = fields[0].to_string();
        let accnum = fields[1].to_string();
        target_lists.entry(accnum).or_default().push(seqid);
    }
    let initial_target_count = target_lists.len();

    // Open accession map
    // Iterate over the accession maps
    // accmap: accnum blank taxid
    let mut accessions_searched = 0;
    for accmap in &args[2..] {
        if target_lists.is_empty() {
            break;
        }

        // Open the current accession map file
        let accmap_file = File::open(accmap)?;
        let reader = BufReader::new(accmap_file);

        // Iterate over the lines in the accession map file
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            let accnum = fields[0].to_string();
            accessions_searched += 1;

            // If the accession number is in the target list, print the sequence ID and taxid
            // Snatch all seqIDs from Accession number in table
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
        // Create lookup file
        let lookup_file_content = "seq1\tacc1\nseq2\tacc2\nseq3\tacc1\n";
        let lookup_file_path = "test_lookup.txt";
        let mut lookup_file = File::create(lookup_file_path).unwrap();
        lookup_file
            .write_all(lookup_file_content.as_bytes())
            .unwrap();

        // Create accmap file
        let accmap_file_content = "acc1\t123\t1000\nacc2\t456\t2000\n";
        let accmap_file_path = "test_accmap.txt";
        let mut accmap_file = File::create(accmap_file_path).unwrap();
        accmap_file
            .write_all(accmap_file_content.as_bytes())
            .unwrap();

        // Call the main function with test arguments
        let _args = vec![
            "lookup_accession_numbers".to_string(),
            lookup_file_path.to_string(),
            accmap_file_path.to_string(),
        ];
        std::env::set_var("RUST_BACKTRACE", "1");
        main().unwrap();
        // Redirect stdout to a file
        let stdout_backup = io::stdout();
        let mut file = File::create("stdout.txt").unwrap();
        {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            handle.write_all(b"Hello, world!").unwrap();
            handle.flush().unwrap();
            file.write_all(&handle.into_inner().unwrap()).unwrap();
        }

        // Now you can read the file and compare its contents to your expected output
        let mut output = String::new();
        file.read_to_string(&mut output).unwrap();
        assert_eq!(output, "seq1\t1000\nseq3\t1000\nseq2\t2000\n");

        // ... your cleanup code ...

        // Clean up test files
        std::fs::remove_file(lookup_file_path).unwrap();
        std::fs::remove_file(accmap_file_path).unwrap();
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
        let _args = vec![
            "lookup_accession_numbers".to_string(),
            lookup_file_path.to_string(),
            accmap_file_path.to_string(),
        ];
        std::env::set_var("RUST_BACKTRACE", "1");
        main().unwrap();

        // Check the output
        let _expected_output = "seq1\t1000\n";
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
