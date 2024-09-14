/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use crate::kraken2_data::{ReadCounter, TaxonCounters, TaxonCounts};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::taxonomy::Taxonomy;
use std::ops::Add;

/// Calculates the clade counts based on the given call counts and taxonomy.
pub fn get_clade_counts(taxonomy: &Taxonomy, call_counts: &TaxonCounts) -> TaxonCounts {
    let mut clade_counts = HashMap::new();

    for (&taxid, &count) in call_counts {
        let mut current_taxid = taxid;
        while current_taxid != 0 {
            *clade_counts.entry(current_taxid).or_insert(0) += count;
            current_taxid =
                taxonomy.nodes[taxonomy.get_internal_id(current_taxid) as usize].parent_id;
        }
    }

    clade_counts
}

/// Calculates the clade counters based on the given call counters and taxonomy.
pub fn get_clade_counters(taxonomy: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let mut clade_counters: HashMap<u64, ReadCounter> = HashMap::new();
    for (&taxid, counter) in call_counters {
        let mut current_taxid = taxid;
        while current_taxid != 0 {
            clade_counters
                .entry(current_taxid)
                .or_default()
                .clone()
                .add(counter.clone()); // Replace .merge(counter) with .add(counter)
            current_taxid =
                taxonomy.nodes[taxonomy.get_internal_id(current_taxid) as usize].parent_id;
        }
    }

    clade_counters
}

/// Writes an MPA-style report line to the given file.
fn write_mpa_style_report_line(
    file: &mut File,
    clade_count: u64,
    taxonomy_line: &str,
) -> std::io::Result<()> {
    writeln!(file, "{}\t{}", taxonomy_line, clade_count)
}

/// Performs a depth-first search to generate the MPA-style report.
fn mpa_report_dfs(
    taxid: u64,
    file: &mut File,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) -> std::io::Result<()> {
    let internal_id = taxonomy.get_internal_id(taxid);
    if !report_zeros && clade_counts.get(&taxid).unwrap_or(&0) == &0 {
        return Ok(());
    }

    let node = &taxonomy.nodes[internal_id as usize];
    let rank = std::str::from_utf8(&taxonomy.rank_data[node.rank_offset as usize..]).unwrap();

    if let Some(rank_code) = get_rank_code(rank) {
        let name = std::str::from_utf8(&taxonomy.name_data[node.name_offset as usize..]).unwrap();
        let name = format!("{}__{}", rank_code, name);
        taxonomy_names.push(name);
        let taxonomy_line = taxonomy_names.join("|");
        write_mpa_style_report_line(
            file,
            clade_counts.get(&taxid).copied().unwrap_or(0),
            &taxonomy_line,
        )?;
    }

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (node.first_child..node.first_child + child_count as u64)
            .map(|internal_id| taxonomy.nodes[internal_id as usize].external_id)
            .collect();
        children.sort_by(|a, b| {
            clade_counts
                .get(b)
                .unwrap_or(&0)
                .cmp(&clade_counts.get(a).unwrap_or(&0))
        });

        for child in children {
            mpa_report_dfs(
                child,
                file,
                report_zeros,
                taxonomy,
                clade_counts,
                taxonomy_names,
            )?;
        }
    }

    if get_rank_code(rank).is_some() {
        taxonomy_names.pop();
    }

    Ok(())
}

/// Generates an MPA-style report and writes it to the specified file.
pub fn report_mpa_style(
    filename: &str,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
) -> std::io::Result<()> {
    let call_counts: TaxonCounts = call_counters
        .iter()
        .map(|(&k, v)| (k, v.read_count()))
        .collect();
    let clade_counts = get_clade_counts(taxonomy, &call_counts);
    let mut file = File::create(filename)?;
    let mut taxonomy_names = Vec::new();
    mpa_report_dfs(
        1,
        &mut file,
        report_zeros,
        taxonomy,
        &clade_counts,
        &mut taxonomy_names,
    )
}

/// Writes a Kraken-style report line to the given file.
fn write_kraken_style_report_line(
    file: &mut File,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &ReadCounter,
    taxon_counter: &ReadCounter,
    rank_str: &str,
    taxid: u64,
    sci_name: &str,
    depth: usize,
) -> std::io::Result<()> {
    let pct = 100.0 * clade_counter.read_count() as f64 / total_seqs as f64;
    write!(
        file,
        "{:6.2}\t{}\t{}\t",
        pct,
        clade_counter.read_count(),
        taxon_counter.read_count()
    )?;

    if report_kmer_data {
        write!(
            file,
            "{}\t{}\t",
            clade_counter.kmer_count(),
            clade_counter.distinct_kmer_count()
        )?;
    }

    write!(file, "{}\t{}\t", rank_str, taxid)?;
    for _ in 0..depth {
        write!(file, "  ")?;
    }
    writeln!(file, "{}", sci_name)
}

/// Performs a depth-first search to generate the Kraken-style report.
fn kraken_report_dfs(
    taxid: u64,
    file: &mut File,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    clade_counters: &TaxonCounters,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    rank_code: char,
    rank_depth: i32,
    depth: usize,
) -> std::io::Result<()> {
    let internal_id = taxonomy.get_internal_id(taxid);
    if !report_zeros
        && clade_counters
            .get(&taxid)
            .map_or(true, |c| c.read_count() == 0)
    {
        return Ok(());
    }

    let node = &taxonomy.nodes[internal_id as usize];
    let rank = std::str::from_utf8(&taxonomy.rank_data[node.rank_offset as usize..]).unwrap();
    let (rank_code, rank_depth) = get_kraken_rank_info(rank, rank_code, rank_depth);
    let rank_str = if rank_depth != 0 {
        format!("{}{}", rank_code, rank_depth)
    } else {
        rank_code.to_string()
    };

    let sci_name = std::str::from_utf8(&taxonomy.name_data[node.name_offset as usize..]).unwrap();
    write_kraken_style_report_line(
        file,
        report_kmer_data,
        total_seqs,
        clade_counters
            .get(&taxid)
            .unwrap_or(&ReadCounter::default()),
        call_counters.get(&taxid).unwrap_or(&ReadCounter::default()),
        &rank_str,
        node.external_id,
        sci_name,
        depth,
    )?;

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (node.first_child..node.first_child + child_count as u64)
            .map(|internal_id| taxonomy.nodes[internal_id as usize].external_id)
            .collect();
        children.sort_by(|a, b| {
            clade_counters
                .get(b)
                .unwrap_or(&ReadCounter::default())
                .read_count()
                .cmp(
                    &clade_counters
                        .get(a)
                        .unwrap_or(&ReadCounter::default())
                        .read_count(),
                )
        });

        for child in children {
            kraken_report_dfs(
                child,
                file,
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

/// Generates a Kraken-style report and writes it to the specified file.
pub fn report_kraken_style(
    filename: &str,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    total_unclassified: u64,
) -> std::io::Result<()> {
    let clade_counters = get_clade_counters(taxonomy, call_counters);
    let mut file = File::create(filename)?;

    if total_unclassified != 0 || report_zeros {
        let rc = ReadCounter::with_counts(total_unclassified, 0);
        write_kraken_style_report_line(
            &mut file,
            report_kmer_data,
            total_seqs,
            &rc,
            &rc,
            "U",
            0,
            "unclassified",
            0,
        )?;
    }

    kraken_report_dfs(
        1,
        &mut file,
        report_zeros,
        report_kmer_data,
        taxonomy,
        &clade_counters,
        call_counters,
        total_seqs,
        'R',
        -1,
        0,
    )
}

/// Retrieves the rank code for the given rank string.
fn get_rank_code(rank: &str) -> Option<char> {
    match rank {
        "superkingdom" => Some('d'),
        "kingdom" => Some('k'),
        "phylum" => Some('p'),
        "class" => Some('c'),
        "order" => Some('o'),
        "family" => Some('f'),
        "genus" => Some('g'),
        "species" => Some('s'),
        _ => None,
    }
}

/// Retrieves the rank information for the Kraken-style report.
fn get_kraken_rank_info(rank: &str, rank_code: char, rank_depth: i32) -> (char, i32) {
    match rank {
        "superkingdom" | "kingdom" | "phylum" | "class" | "order" | "family" | "genus"
        | "species" => (get_rank_code(rank).unwrap(), 0),
        _ => (rank_code, rank_depth + 1),
    }
}

#[cfg(test)]
mod tests {
    use crate::hyperloglogplus::murmurhash3_finalizer;

    use super::*;

    #[test]
    fn test_get_clade_counts() {
        let mut taxonomy = Taxonomy::new();
        taxonomy.nodes = vec![
            TaxonomyNode {
                parent_id: 0,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 1,
                godparent_id: 0,
            },
            TaxonomyNode {
                parent_id: 1,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 2,
                godparent_id: 0,
            },
            TaxonomyNode {
                parent_id: 2,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 3,
                godparent_id: 0,
            },
        ];
        taxonomy.generate_external_to_internal_id_map();

        let call_counts = vec![(1, 10), (2, 5), (3, 3)].into_iter().collect();
        let clade_counts = get_clade_counts(&taxonomy, &call_counts);

        assert_eq!(clade_counts.get(&1), Some(&18)); // I get Some(&10)
        assert_eq!(clade_counts.get(&2), Some(&8));
        assert_eq!(clade_counts.get(&3), Some(&3));
    }

    #[test]
    fn test_get_clade_counters() {
        let mut taxonomy = Taxonomy::new();
        taxonomy.nodes = vec![
            TaxonomyNode {
                parent_id: 0,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 1,
                godparent_id: 0,
            },
            TaxonomyNode {
                parent_id: 1,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 2,
                godparent_id: 0,
            },
            TaxonomyNode {
                parent_id: 2,
                first_child: 0,
                child_count: 0,
                name_offset: 0,
                rank_offset: 0,
                external_id: 3,
                godparent_id: 0,
            },
        ];
        taxonomy.generate_external_to_internal_id_map();

        let call_counters = vec![
            (
                1,
                ReadCounter::with_counts(
                    10,
                    100,
                    HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer),
                ),
            ),
            (
                2,
                ReadCounter::with_counts(
                    5,
                    50,
                    HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer),
                ),
            ),
            (
                3,
                ReadCounter::with_counts(
                    3,
                    30,
                    HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer),
                ),
            ),
        ]
        .into_iter()
        .collect();
        let clade_counters = get_clade_counters(&taxonomy, &call_counters);

        assert_eq!(clade_counters.get(&1).unwrap().read_count(), 18);
        assert_eq!(clade_counters.get(&2).unwrap().read_count(), 8);
        assert_eq!(clade_counters.get(&3).unwrap().read_count(), 3);
    }

    #[test]
    fn test_get_rank_code() {
        assert_eq!(get_rank_code("superkingdom"), Some('d'));
        assert_eq!(get_rank_code("kingdom"), Some('k'));
        assert_eq!(get_rank_code("phylum"), Some('p'));
        assert_eq!(get_rank_code("class"), Some('c'));
        assert_eq!(get_rank_code("order"), Some('o'));
        assert_eq!(get_rank_code("family"), Some('f'));
        assert_eq!(get_rank_code("genus"), Some('g'));
        assert_eq!(get_rank_code("species"), Some('s'));
        assert_eq!(get_rank_code("subspecies"), None);
    }
}
