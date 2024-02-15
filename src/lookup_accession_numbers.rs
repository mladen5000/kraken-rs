use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lookup_accession_numbers <lookup file> <accmaps>");
        std::process::exit(1);
    }

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

    for accmap in &args[2..] {
        if target_lists.is_empty() {
            break;
        }

        let accmap_file = File::open(accmap)?;
        let reader = BufReader::new(accmap_file);

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
