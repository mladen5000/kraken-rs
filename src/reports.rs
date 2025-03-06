use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::str;

use crate::kraken2_data::{TaxId, TaxonCounters, TaxonCountersMap, TaxonCounts};
use crate::taxonomy::Taxonomy;

pub fn get_clade_counts(tax: &Taxonomy, call_counts: &TaxonCounts) -> TaxonCounts {
    let mut clade_counts = TaxonCounts::new();
    for (&taxid, &count) in call_counts {
        let mut current = taxid;
        while current != 0 {
            *clade_counts.entry(current).or_insert(0) += count;
            current = tax.nodes()[current as usize].parent_id;
        }
    }
    clade_counts
}

pub fn get_clade_counters(tax: &Taxonomy, call_counters: &TaxonCountersMap) -> TaxonCountersMap {
    let mut clade_counters = TaxonCountersMap::new();
    for (&taxid, counter) in call_counters {
        let mut current = taxid;
        while current != 0 {
            let entry = clade_counters
                .entry(current)
                .or_insert_with(|| TaxonCounters::new_with_precision(10));
            entry.merge(counter);
            current = tax.nodes()[current as usize].parent_id;
        }
    }
    clade_counters
}

fn print_mpa_style_report_line(
    ofs: &mut BufWriter<File>,
    clade_count: u64,
    taxonomy_line: &str,
) -> Result<()> {
    writeln!(ofs, "{}\t{}", taxonomy_line, clade_count)?;
    Ok(())
}

fn mpa_report_dfs(
    taxid: TaxId,
    ofs: &mut BufWriter<File>,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) -> Result<()> {
    if !report_zeros && clade_counts.get(&taxid).copied().unwrap_or(0) == 0 {
        return Ok(());
    }

    let node = &taxonomy.nodes()[taxid as usize];
    let rank = taxonomy.rank_data_at(node.rank_offset);

    let rank_code = match rank {
        "superkingdom" => 'd',
        "kingdom" => 'k',
        "phylum" => 'p',
        "class" => 'c',
        "order" => 'o',
        "family" => 'f',
        "genus" => 'g',
        "species" => 's',
        _ => '\0',
    };

    if rank_code != '\0' {
        let name = taxonomy.name_data_at(node.name_offset);
        let tax_name = format!("{}__{}", rank_code, name);
        taxonomy_names.push(tax_name);

        let taxonomy_line = taxonomy_names.join("|");
        print_mpa_style_report_line(ofs, clade_counts[&taxid], &taxonomy_line)?;
    }

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (0..child_count)
            .map(|i| node.first_child + i as u64)
            .collect();
        children.sort_by(|&a, &b| {
            let ca = clade_counts.get(&a).copied().unwrap_or(0);
            let cb = clade_counts.get(&b).copied().unwrap_or(0);
            cb.cmp(&ca)
        });

        for child in children {
            mpa_report_dfs(
                child,
                ofs,
                report_zeros,
                taxonomy,
                clade_counts,
                taxonomy_names,
            )?;
        }
    }

    if rank_code != '\0' {
        taxonomy_names.pop();
    }

    Ok(())
}

pub fn report_mpa_style(
    filename: &str,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCountersMap,
) -> Result<()> {
    let mut call_counts = TaxonCounts::new();
    for (&taxid, counter) in call_counters {
        call_counts.insert(taxid, counter.get_read_count());
    }
    let clade_counts = get_clade_counts(taxonomy, &call_counts);

    let file = File::create(filename).with_context(|| format!("error creating {}", filename))?;
    let mut ofs = BufWriter::new(file);

    let mut taxonomy_names = Vec::new();
    mpa_report_dfs(
        1,
        &mut ofs,
        report_zeros,
        taxonomy,
        &clade_counts,
        &mut taxonomy_names,
    )?;
    Ok(())
}

fn print_kraken_style_report_line(
    ofs: &mut BufWriter<File>,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &TaxonCounters,
    taxon_counter: &TaxonCounters,
    rank_str: &str,
    taxid: u32,
    sci_name: &str,
    depth: i32,
) -> Result<()> {
    write!(
        ofs,
        "{:.2}\t{}\t{}",
        100.0 * (clade_counter.get_read_count() as f64) / (total_seqs as f64),
        clade_counter.get_read_count(),
        taxon_counter.get_read_count()
    )?;

    if report_kmer_data {
        write!(
            ofs,
            "\t{}\t{}",
            clade_counter.get_read_count(),
            clade_counter.get_kmer_distinct()
        )?;
    }

    write!(ofs, "\t{}\t{}\t", rank_str, taxid)?;
    for _ in 0..depth {
        write!(ofs, "  ")?;
    }
    writeln!(ofs, "{}", sci_name)?;

    Ok(())
}

fn kraken_report_dfs(
    taxid: u32,
    ofs: &mut BufWriter<File>,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    clade_counters: &TaxonCountersMap,
    call_counters: &TaxonCountersMap,
    total_seqs: u64,
    rank_code: char,
    rank_depth: i32,
    depth: i32,
) -> Result<()> {
    let taxid_u64 = taxid as u64;
    if !report_zeros
        && clade_counters
            .get(&taxid_u64)
            .map_or(0, |c| c.get_read_count())
            == 0
    {
        return Ok(());
    }

    let node = &taxonomy.nodes()[taxid as usize];
    let rank = taxonomy.rank_data_at(node.rank_offset);

    let (rank_code, rank_depth) = match rank {
        "superkingdom" => ('D', 0),
        "kingdom" => ('K', 0),
        "phylum" => ('P', 0),
        "class" => ('C', 0),
        "order" => ('O', 0),
        "family" => ('F', 0),
        "genus" => ('G', 0),
        "species" => ('S', 0),
        _ => (rank_code, rank_depth + 1),
    };

    let rank_str = if rank_depth != 0 {
        format!("{}{}", rank_code, rank_depth)
    } else {
        rank_code.to_string()
    };

    let sci_name = taxonomy.name_data_at(node.name_offset);
    let empty_counter = TaxonCounters::new_with_precision(10);
    let taxon_counter = call_counters.get(&taxid_u64).unwrap_or(&empty_counter);
    let clade_counter = clade_counters.get(&taxid_u64).unwrap_or(&empty_counter);

    print_kraken_style_report_line(
        ofs,
        report_kmer_data,
        total_seqs,
        clade_counter,
        taxon_counter,
        &rank_str,
        node.external_id as u32,
        sci_name,
        depth,
    )?;

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (0..child_count)
            .map(|i| node.first_child + i as u64)
            .collect();
        children.sort_by(|&a, &b| {
            let ca = clade_counters.get(&a).map_or(0, |c| c.get_read_count());
            let cb = clade_counters.get(&b).map_or(0, |c| c.get_read_count());
            cb.cmp(&ca)
        });

        for child in children {
            kraken_report_dfs(
                child as u32,
                ofs,
                report_zeros,
                report_kmer_data,
                taxonomy,
                clade_counters,
                call_counters,
                total_seqs,
                rank_code,
                rank_depth,
                depth + 1,
            )?;
        }
    }

    Ok(())
}

pub fn report_kraken_style(
    filename: &str,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCountersMap,
    total_seqs: u64,
    total_unclassified: u64,
) -> Result<()> {
    let file = File::create(filename).with_context(|| format!("Error creating {}", filename))?;
    let mut ofs = BufWriter::new(file);

    // Handle unclassified sequences
    if total_unclassified != 0 || report_zeros {
        let mut unclass_counter = TaxonCounters::new_with_precision(10);
        for _ in 0..total_unclassified {
            unclass_counter.increment_read_count();
        }
        print_kraken_style_report_line(
            &mut ofs,
            report_kmer_data,
            total_seqs,
            &unclass_counter,
            &unclass_counter,
            "U",
            0,
            "unclassified",
            0,
        )?;
    }

    let clade_counters = get_clade_counters(taxonomy, call_counters);

    // DFS from root (taxid=1)
    kraken_report_dfs(
        1,
        &mut ofs,
        report_zeros,
        report_kmer_data,
        taxonomy,
        &clade_counters,
        call_counters,
        total_seqs,
        'R',
        -1,
        0,
    )?;

    Ok(())
}

pub fn report_stats(
    duration_secs: f64,
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
) {
    if atty::is(atty::Stream::Stderr) {
        eprint!("\r");
    }

    let total_unclassified = total_sequences - total_classified;

    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({} Kseq/m, {:.2} Mbp/m).",
        total_sequences,
        total_bases as f64 / 1.0e6,
        duration_secs,
        total_sequences as f64 / 1.0e3 / (duration_secs / 60.0),
        total_bases as f64 / 1.0e6 / (duration_secs / 60.0)
    );

    eprintln!(
        "  {} sequences classified ({:.2}%)",
        total_classified,
        total_classified as f64 * 100.0 / total_sequences as f64
    );

    eprintln!(
        "  {} sequences unclassified ({:.2}%)",
        total_unclassified,
        total_unclassified as f64 * 100.0 / total_sequences as f64
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_report_kraken_style() {
        // Create a mock taxonomy and counts
        let taxonomy = Taxonomy::default();
        let mut counts = TaxonCountersMap::new();

        // Add some test data
        let counter = TaxonCounters::new_with_precision(10);
        counts.insert(1, counter);

        // Test report generation
        let result =
            report_kraken_style("test_kraken.txt", true, false, &taxonomy, &counts, 100, 10);
        assert!(result.is_ok());
    }

    #[test]
    fn test_report_mpa_style() {
        // Create a mock taxonomy and counts
        let taxonomy = Taxonomy::default();
        let mut counts = TaxonCountersMap::new();

        // Add some test data
        let counter = TaxonCounters::new_with_precision(10);
        counts.insert(1, counter);

        // Test report generation
        let result = report_mpa_style("test_mpa.txt", true, &taxonomy, &counts);
        assert!(result.is_ok());
    }
}
