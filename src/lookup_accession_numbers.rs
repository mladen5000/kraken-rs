use std::collections::HashMap;
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::process::exit;

fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lookup_accession_numbers <lookup file> <accmaps>");
        exit(1);
    }

    let lookup_file = &args[1];
    let accmap_files = &args[2..];

    // Create a map to store target accession numbers and their associated sequence IDs
    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();

    // Open and read the lookup file
    let lookup_list_file = File::open(lookup_file).expect("Failed to open lookup file");
    let reader = BufReader::new(lookup_list_file);

    // Parse the lookup file and populate the target_lists map
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let fields: Vec<&str> = line.splitn(2, '\t').collect();
        if fields.len() < 2 {
            continue; // Skip lines that don't have at least two fields
        }
        let seqid = fields[0].to_string();
        let accnum = fields[1].to_string();
        target_lists.entry(accnum).or_default().push(seqid);
    }

    let initial_target_count = target_lists.len();
    let mut accessions_searched: u64 = 0;

    // Display initial progress if stderr is a TTY
    if atty::is(atty::Stream::Stderr) {
        eprint!("\rFound 0/{} targets...", initial_target_count);
    }

    // Iterate over the accession map files provided as arguments
    for accmap_filename in accmap_files {
        if target_lists.is_empty() {
            break; // Stop if all targets have been found
        }

        let accmap_file = File::open(accmap_filename).expect("Failed to open accmap file");
        let reader = BufReader::new(accmap_file);
        let mut lines = reader.lines();

        // Skip header line
        lines.next();

        // Process each line in the accession map file
        for line_result in lines {
            let line = line_result.expect("Failed to read line");
            accessions_searched += 1;

            let mut fields_iter = line.split('\t');
            let accnum = match fields_iter.next() {
                Some(value) => value.to_string(),
                None => {
                    eprintln!("Expected TAB not found in {}", accmap_filename);
                    break;
                }
            };

            // Check if the accession number is in our target list
            if target_lists.contains_key(&accnum) {
                // Skip the next two fields to reach the taxid
                let taxid = fields_iter.nth(2).unwrap_or_default().to_string();

                // Output the sequence IDs and taxid
                if let Some(seqids) = target_lists.remove(&accnum) {
                    for seqid in seqids {
                        println!("{}\t{}", seqid, taxid);
                    }
                }

                // Update progress display
                if atty::is(atty::Stream::Stderr) {
                    eprint!(
                        "\rFound {}/{} targets, searched through {} accession IDs...",
                        initial_target_count - target_lists.len(),
                        initial_target_count,
                        accessions_searched
                    );
                }

                if target_lists.is_empty() {
                    break; // All targets found
                }
            }

            // Periodically update progress display
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

    // Clear the progress display
    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }

    // Final progress report
    eprintln!(
        "Found {}/{} targets, searched through {} accession IDs, search complete.",
        initial_target_count - target_lists.len(),
        initial_target_count,
        accessions_searched
    );

    // Handle any accession numbers that remain unmapped
    if !target_lists.is_empty() {
        eprintln!(
            "lookup_accession_numbers: {}/{} accession numbers remain unmapped, see unmapped.txt",
            target_lists.len(),
            initial_target_count
        );
        let mut ofs = File::create("unmapped.txt").expect("Failed to create unmapped.txt");
        for accnum in target_lists.keys() {
            writeln!(ofs, "{}", accnum).expect("Failed to write to unmapped.txt");
        }
    }
}
