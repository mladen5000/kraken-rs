use std::{
    collections::{HashMap, HashSet}, io::Write, rc
};

// Enumerations
use crate::ncbi_taxonomy;
enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}

type SparseListType = HashSet<u32>;

// Structures
struct TaxonomyNode {
    parent_id: u32,
    rank_offset: usize,
    name_offset: usize,
    external_id: u32,
    child_count: u32,
    first_child: u32,
    // Additional fields as needed
}

struct Taxonomy {
    // Fields and methods related to taxonomy data
}
struct HyperLogLogPlusMinus<T> {
    p: u8,
    m: u8,
    M: Vec<u8>,
    n_observed: u64, // your code here
    sparse: bool,
    sparse_list: SparseListType,
    bit_mixera: u64,
}

struct ReadCounts<T> {
    n_reads: u64,
    n_kmers: u64,
    kmers: T, // distinct kmer count per taxon,
}

type ReadCounter = ReadCounts<HyperLogLogPlusMinus<u64>>;

struct TaxonCounts(HashMap<u32, u64>); // Map of taxon ID to count
struct TaxonCounters(HashMap<u32, ReadCounter>); // Map of taxon ID to ReadCounter
                                                 // type TaxonCounters = HashMap<u32, ReadCounter>; // Map of taxon ID to ReadCounter

// Functions
fn get_clade_counts(tax: &Taxonomy, call_counts: &TaxonCounts) -> TaxonCounts {
    let clade_counts: TaxonCounts = TaxonCounts(HashMap::new());

    for (taxid, count) in call_counts.0 {
        while taxid != 0 {
            let val: &mut u64 = clade_counts.0.entry(taxid).or_insert(0);
            *val += count;
            taxid = tax.nodes()[taxid].parent_id;
        }
    }
    clade_counts
}

fn get_clade_counters(tax: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let clade_counters: TaxonCounters;

    for (taxid, count) in call_counters.0 {
        while taxid != 0 {
            let val = clade_counters.0.entry(taxid).or_insert(ReadCounter {0,0,""});
            // *val += count;
            taxid = tax.nodes()[taxid].parent_id;
        }
    }
    clade_counters
}

fn print_mpa_style_report_line<W: std::io::Write>(
    writer: &mut W,
    clade_count: u64,
    taxonomy_line: &str,
) {
    write!(writer, "{}\t{}", taxonomy_line, clade_count);
}

fn mpa_report_dfs<W: Write>(
    taxid: u32,
    writer: &mut W,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) {
    if !report_zeros && clade_counts[taxid] == 0 {
        return;
    }
    let node = taxonomy.nodes()[taxid];
    let rank = taxonomy.rank_data() + node.rank_offset;

    let rank_code: char = '\0';
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
        let mut name = String::new();
        name.push(rank_code);
        name.push_str("__");
        name.push_str(taxonomy.name_data() + node.name_offset);
        taxonomy_names.push(name);
        let taxonomy_line = "";
        for str in taxonomy_names {
            taxonomy_line.to_string().push_str(&format!("{}|", str));
        }
        taxonomy_line.pop();
        print_mpa_style_report_line(&mut writer, clade_counts[taxid], taxonomy_line);
    }
    let child_count = node.child_count;
    let mut children: Vec<u64>;
    if child_count != 0 {
        let mut children = child_count;
    }
    for i in (0..child_count) {
        children.insert(i, node.first_child + i)
    }
    children.sort_by(|&a, &b| clade_counts[&b].cmp(&clade_counts[&a]));

    for child in children {
        mpa_report_dfs(
            child,
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

fn report_mpa_style<W: Write>(
    writer: &mut W,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
) {
    let mut call_counts: TaxonCounts = HashMap::new();
    for kv_pair in call_counters {
        call_counts[kv_pair.0] = kv_pair.1.read_count()
    }
    let clade_counts = get_clade_counts(taxonomy, &call_counts);
    // Implementation
}

fn print_kraken_style_report_line<W: Write>(
    writer: &mut W,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &ReadCounter,
    taxon_counter: &ReadCounter,
    rank_str: &str,
    taxid: u32,
    sci_name: &str,
    depth: usize,
) {
    let total_seqs: u64;
    let (clade_counter, taxon_counter): (ReadCounter, ReadCounter);
    let depth: i64;
    let taxid: u32;
    let rank_str = String::new();
    let sci_name = String::new();
    let pct_buffer = format!(
        "{{:6.2}}",
        100.0 * clade_counter.read_count() as f64 / total_seqs as f64
    );
    write!(
        writer,
        "{}\\t{}\\t{}\\t",
        pct_buffer,
        clade_counter.read_count(),
        taxon_counter.read_count()
    );
    if report_kmer_data {
        write!(writer,"{clade_counter.kmer_count()}\t{clade_counter.distinct_kmer_count()}\t");
    }
    write!(writer, "{rank_str}\t{taxid}\t");
    for i in (0..depth) {
        write!(writer," ");
    }
    write!(sci_name);

fn kraken_report_dfs<W: Write>(
    taxid: u32,
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
    let clade_counters: TaxonCounters;
    let call_counters: TaxonCounters;
    let total_seqs: u64;
    let rank_code: char;
    let (rank_depth, depth): (i64,i64) = (0,0);

    if ! report_zeros && clade_counters[taxid].read_count() == 0 { return;}
    let node: TaxonomyNode = taxonomy.nodes()[taxid];

    match taxonomy.rank_data() + node.rank_offset {
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
        _ => rank_depth+=1
    }
    // Implementation
}

rank_str(&rank_code, 0, 1);
if rank_depth != 0 {
    rank_str += rank_depth.to_string()
}
let name = taxonomy.name_data() + node.name_offset;

fn report_kraken_style<W: Write>(
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
    let (total_seqs, total_unclassified): (u64,u64) = (0, 0);

    let clade_counters: TaxonCounters = get_clade_counters(&taxonomy, &call_counters);

    let mut ofs = File::create(filename).expect("Unable to create file");

    if total_unclassified != 0 || report_zeros {
        print_kraken_style_report_line(writer, report_kmer_data, total_seqs, rc, rc, "U", 0, "unclassified", )
    }

    kraken_report_dfs(taxid, writer, report_zeros, report_kmer_data, &taxonomy, &clade_counters, &call_counters, total_seqs, 'R', -1, 0)

}

// Unit Tests
#[cfg(test)]
mod tests {
    // Test implementations
}
