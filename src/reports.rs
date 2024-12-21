use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write}; // or READCOUNTER if differently named

use crate::kraken2_data::taxid_t;
use crate::readcounts::TaxonCount;
use crate::taxonomy::{Taxonomy, TaxonomyNode};

type TaxonCounts = HashMap<taxid_t, u32>;
type TaxonCounters = HashMap<taxid_t, TaxonCount>;

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

pub fn get_clade_counters(tax: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let mut clade_counters = TaxonCounters::new();
    for (&taxid, counter) in call_counters {
        let mut current = taxid;
        while current != 0 {
            clade_counters.entry(current).or_default().add(counter);
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
    taxid: taxid_t,
    ofs: &mut BufWriter<File>,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) -> Result<()> {
    let clade_count = clade_counts.get(&taxid).copied().unwrap_or(0);
    if !report_zeros && clade_count == 0 {
        return Ok(());
    }

    let node = taxonomy.nodes()[taxid as usize];
    let rank = taxonomy.rank_data_at(node.rank_offset);

    let rank_code = match rank.as_str() {
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

    let mut added = false;
    if rank_code != '\0' {
        let name = taxonomy.name_data_at(node.name_offset);
        let mut tax_name = String::new();
        tax_name.push(rank_code);
        tax_name.push_str("__");
        tax_name.push_str(name);
        taxonomy_names.push(tax_name);

        let taxonomy_line = taxonomy_names.join("|");
        print_mpa_style_report_line(ofs, clade_count as u64, &taxonomy_line)?;
        added = true;
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

    if added {
        taxonomy_names.pop();
    }

    Ok(())
}

pub fn report_mpa_style(
    filename: &str,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
) -> Result<()> {
    let mut call_counts = TaxonCounts::new();
    for (&taxid, counter) in call_counters {
        call_counts.insert(taxid, counter.read_count());
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
    clade_counter: &TaxonCount,
    taxon_counter: &TaxonCount,
    rank_str: &str,
    taxid: u32,
    sci_name: &str,
    depth: i32,
) -> Result<()> {
    let pct = 100.0 * (clade_counter.read_count() as f64) / (total_seqs as f64);

    write!(
        ofs,
        "{:6.2}\t{}\t{}\t",
        pct,
        clade_counter.read_count(),
        taxon_counter.read_count()
    )?;
    if report_kmer_data {
        write!(
            ofs,
            "{}\t{}\t",
            clade_counter.kmer_count(),
            clade_counter.distinct_kmer_count()
        )?;
    }

    write!(ofs, "{}\t{}\t", rank_str, taxid)?;
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
    clade_counters: &TaxonCounters,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    mut rank_code: char,
    mut rank_depth: i32,
    depth: i32,
) -> Result<()> {
    let clade_counter = clade_counters
        .get(&(taxid as taxid_t))
        .unwrap_or(&TaxonCount::default());
    if !report_zeros && clade_counter.read_count() == 0 && taxid != 1 {
        return Ok(());
    }

    let node = taxonomy.nodes()[taxid as usize];
    let rank = taxonomy.rank_data_at(node.rank_offset);

    match rank.as_str() {
        "superkingdom" => {
            rank_code = 'D';
            rank_depth = 0;
        }
        "kingdom" => {
            rank_code = 'K';
            rank_depth = 0;
        }
        "phylum" => {
            rank_code = 'P';
            rank_depth = 0;
        }
        "class" => {
            rank_code = 'C';
            rank_depth = 0;
        }
        "order" => {
            rank_code = 'O';
            rank_depth = 0;
        }
        "family" => {
            rank_code = 'F';
            rank_depth = 0;
        }
        "genus" => {
            rank_code = 'G';
            rank_depth = 0;
        }
        "species" => {
            rank_code = 'S';
            rank_depth = 0;
        }
        _ => {
            rank_depth += 1;
        }
    }

    let mut rank_str = rank_code.to_string();
    if rank_depth != 0 {
        rank_str.push_str(&rank_depth.to_string());
    }

    let sci_name = taxonomy.name_data_at(node.name_offset);
    let taxon_counter = call_counters
        .get(&(taxid as taxid_t))
        .unwrap_or(&TaxonCount::default());

    print_kraken_style_report_line(
        ofs,
        report_kmer_data,
        total_seqs,
        clade_counter,
        taxon_counter,
        &rank_str,
        node.external_id,
        sci_name,
        depth,
    )?;

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (0..child_count)
            .map(|i| node.first_child + i as u64)
            .collect();
        children.sort_by(|&a, &b| {
            let ca = clade_counters.get(&a).map(|x| x.read_count()).unwrap_or(0);
            let cb = clade_counters.get(&b).map(|x| x.read_count()).unwrap_or(0);
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
    call_counters: &TaxonCounters,
    total_seqs: u64,
    total_unclassified: u64,
) -> Result<()> {
    let clade_counters = get_clade_counters(taxonomy, call_counters);

    let file = File::create(filename).with_context(|| format!("Error creating {}", filename))?;
    let mut ofs = BufWriter::new(file);

    // Handle unclassified
    if total_unclassified != 0 || report_zeros {
        let mut rc = TaxonCount::default();
        rc.set_read_count(total_unclassified);
        print_kraken_style_report_line(
            &mut ofs,
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
