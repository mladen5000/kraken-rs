use std::{
    collections::{HashMap, HashSet},
    io::Write,
};

use crate::taxonomy::{Taxonomy, TaxonomyNode};
// Enumerations
enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}

type SparseListType = HashSet<u64>;

pub struct HyperLogLogPlusMinus {
    p: u8,
    m: u8,
    M: Vec<u8>,
    n_observed: u64, // your code here
    sparse: bool,
    sparse_list: SparseListType,
    bit_mixera: u64,
}

pub struct ReadCounts<T> {
    n_reads: u64,
    n_kmers: u64,
    kmers: T, // distinct kmer count per taxon,
}
impl<T> ReadCounts<T> {
    fn new(n_reads: u64, n_kmers: u64, kmers: T) -> Self {
        Self {
            n_reads,
            n_kmers,
            kmers,
        }
    }
    fn increment_read_counts(&mut self) {
        self.n_reads += 1;
    }
    fn kmer_counts(&mut self) -> u64 {
        self.n_kmers
    }
    fn distinct_kmer_counts(&mut self) {
        unimplemented!()
    }
}

impl<T> ReadCounts<T> {
    pub fn read_count(&self) -> u64 {
        self.n_reads
    }
    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }
    pub fn kmer_counter(&self) -> &T {
        &self.kmers
    }
    pub fn distinct_kmer_count(&self) -> u64 {
        todo!()
    }
}

type ReadCounter = ReadCounts<HyperLogLogPlusMinus>;

pub struct TaxonCounts(pub HashMap<u64, u64>); // Map of taxon ID to count
impl TaxonCounts {
    pub fn new() -> Self {
        TaxonCounts { 0: HashMap::new() }
    }
    pub fn iter(&self) -> std::collections::hash_map::Iter<u64, u64> {
        self.0.iter()
    }

    pub fn iter_mut(&mut self) -> std::collections::hash_map::IterMut<u64, u64> {
        self.0.iter_mut()
    }
}
impl Iterator for TaxonCounts {
    type Item = (u64, u64);
    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}
pub struct TaxonCounters(HashMap<u64, ReadCounter>); // Map of taxon ID to ReadCounter

impl TaxonCounters {
    pub fn new() -> Self {
        TaxonCounters { 0: HashMap::new() }
    }
}
impl Default for TaxonCounters {
    fn default() -> Self {
        Self::new()
    }
}
// type TaxonCounters = HashMap<u64, ReadCounter>; // Map of taxon ID to ReadCounter

// Functions
fn get_clade_counts(tax: &Taxonomy, call_counts: &TaxonCounts) -> TaxonCounts {
    let mut clade_counts: TaxonCounts = TaxonCounts(HashMap::new());

    for (mut taxid, count) in call_counts.0 {
        while taxid != 0 {
            let val: &mut u64 = clade_counts.0.entry(taxid).or_insert(0);
            *val += count;
            taxid = tax.nodes()[taxid as usize].parent_id;
        }
    }
    clade_counts
}

fn get_clade_counters(tax: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let mut clade_counters: TaxonCounters;

    for (mut taxid, count) in call_counters.0 {
        while taxid != 0 {
            let val = clade_counters.0.entry(taxid).or_insert(ReadCounter {
                n_reads: 0,
                n_kmers: 0,
                kmers: HyperLogLogPlusMinus {
                    p: 0,
                    m: 0,
                    M: Vec::new(),
                    n_observed: 0,
                    sparse: false,
                    sparse_list: HashSet::new(),
                    bit_mixera: 0,
                },
            });
            // *val += count;
            taxid = tax.nodes()[taxid as usize].parent_id;
        }
    }
    clade_counters
}

pub fn print_mpa_style_report_line<W: std::io::Write>(
    writer: &mut W,
    clade_count: u64,
    taxonomy_line: &str,
) {
    write!(writer, "{}\t{}", taxonomy_line, clade_count);
}

pub fn mpa_report_dfs<W: Write>(
    taxid: usize,
    writer: &mut W,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) {
    let taxref = taxid as u64;
    if !report_zeros && clade_counts.0[&taxref] == 0 {
        return;
    }
    let node = taxonomy.nodes()[taxid];
    let rank = format!("{}{}", taxonomy.rank_data(), &node.rank_offset).as_ref();

    let mut rank_code: char = '\0';
    rank_code = match rank {
        "superkingdom" => 'd',
        "kingdom" => 'k',
        "phylum" => 'p',
        "class" => 'c',
        "order" => 'o',
        "family" => 'f',
        "genus" => 'g',
        "species" => 's',
    };

    if rank_code.is_alphanumeric() {
        let name = format!(
            "{}__{}{}",
            rank_code,
            taxonomy.name_data(),
            node.name_offset
        );
        taxonomy_names.push(name);
        let mut taxonomy_line = String::new();
        for str in taxonomy_names {
            taxonomy_line.to_string().push_str(&format!("{}|", str));
        }
        taxonomy_line.pop();
        let taxref = taxid as u64;
        print_mpa_style_report_line(writer, clade_counts.0[&taxref], &taxonomy_line);
    }
    let child_count = node.child_count;
    let mut children: Vec<u64>;
    if child_count != 0 {
        let mut children = child_count;
    }
    for i in 0..child_count {
        children.insert(i as usize, node.first_child + i)
    }
    children.sort_by(|&a, &b| clade_counts.0[&b].cmp(&clade_counts.0[&a]));

    for child in children {
        mpa_report_dfs(
            child as usize,
            writer,
            report_zeros,
            taxonomy,
            clade_counts, // Dereference the &TaxonCounts to access the underlying TaxonCounts value
            taxonomy_names,
        );
    }

    if rank_code.is_alphanumeric() {
        taxonomy_names.pop();
    }

    // Implementation
}

pub fn report_mpa_style<W: Write>(
    writer: &mut W,
    report_zeros: bool,
    taxonomy: &mut Taxonomy,
    call_counters: &TaxonCounters,
) {
    let mut call_counts = TaxonCounts(HashMap::new());
    for kv_pair in call_counters.0 {
        call_counts.0[&kv_pair.0] = kv_pair.1.read_count()
    }
    let clade_counts = get_clade_counts(taxonomy, &call_counts);
    // Implementation
}

pub fn print_kraken_style_report_line<W: Write>(
    writer: &mut W,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &ReadCounter,
    taxon_counter: &ReadCounter,
    rank_str: &str,
    taxid: u64,
    sci_name: &str,
    depth: usize,
) {
    let total_seqs: u64;
    let (clade_counter, taxon_counter): (ReadCounter, ReadCounter);
    let depth: usize;
    let taxid: u64;
    let rank_str = String::new();
    let sci_name = String::new();

    let percent_value = 100.0 * clade_counter.read_count() as f64 / total_seqs as f64;
    let pct_buffer = format!("{{percent_value:6.2}}",);
    write!(
        writer,
        "{}\\t{}\\t{}\\t",
        pct_buffer,
        clade_counter.read_count(),
        taxon_counter.read_count()
    );
    if report_kmer_data {
        write!(
            writer,
            "{}\t{}\t",
            clade_counter.kmer_count(),
            clade_counter.distinct_kmer_count()
        );
    }
    write!(writer, "{rank_str}\t{taxid}\t");
    for i in 0..depth {
        write!(writer, " ");
    }
    write!(writer, "{sci_name}");
}

fn kraken_report_dfs<W: Write>(
    taxid: u64,
    writer: &mut W,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    clade_counters: &TaxonCounters,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    rank_code: char,
    rank_depth: i32,
    depth: usize,
) {
    let report_kmer_data: bool;
    let taxonomy: Taxonomy;
    let clade_counters = TaxonCounters(HashMap::new());
    let call_counters = TaxonCounters(HashMap::new());
    let total_seqs: u64;
    let rank_code: char;
    let (mut rank_depth, depth): (usize, usize) = (0, 0);

    if !report_zeros && clade_counters.0[&taxid].read_count() == 0 {
        return;
    }
    let node: TaxonomyNode = taxonomy.nodes()[taxid as usize];

    match format!("{}{}", taxonomy.rank_data(), node.rank_offset).as_ref() {
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
        _ => rank_depth += 1,
    }
    // Implementation

    let mut rank_str = format!("{}", &rank_code);
    if rank_depth != 0 {
        rank_str.push_str(&rank_depth.to_string());
    }
    let name = format!("{}{}", taxonomy.name_data(), node.name_offset);

    print_kraken_style_report_line(
        writer,
        report_kmer_data,
        total_seqs,
        &clade_counters.0[&taxid],
        &call_counters.0[&taxid],
        &rank_str,
        taxid,
        &name,
        depth,
    );
}

pub fn report_kraken_style<W: Write>(
    writer: &mut W,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    total_unclassified: u64,
) {
    let taxonomy: Taxonomy;
    let call_counters: TaxonCounters;
    let (total_seqs, total_unclassified): (u64, u64) = (0, 0);

    let clade_counters: TaxonCounters = get_clade_counters(&taxonomy, &call_counters);

    if total_unclassified != 0 || report_zeros {
        let rc: ReadCounts<HyperLogLogPlusMinus> = ReadCounter::new(
            total_unclassified,
            0,
            HyperLogLogPlusMinus {
                p: 0,
                m: 0,
                M: Vec::new(),
                n_observed: 0,
                sparse: false,
                sparse_list: HashSet::new(),
                bit_mixera: 0,
            },
        );
        print_kraken_style_report_line(
            writer,
            report_kmer_data,
            total_seqs,
            &rc,
            &rc,
            "U",
            0,
            "unclassified",
            0,
        )
    }

    kraken_report_dfs(
        1,
        writer,
        report_zeros,
        report_kmer_data,
        &taxonomy,
        &clade_counters,
        &call_counters,
        total_seqs,
        'R',
        -1,
        0,
    )
}

// Additional reporting functions as needed
// Unit Tests
#[cfg(test)]
mod tests {
    // Test implementations
}
