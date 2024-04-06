// pub mod aa_translate;
// pub mod classify;
// pub mod classify3;
// pub mod estimate_capacity;
pub mod hyperloglogplus;

pub mod kraken2_data;
pub mod kv_store;
// pub mod readcounts_old;
pub mod reports;
pub mod threadpool;

// Build DB
pub mod build_db;
// pub mod compact_hash;
pub mod mmap_file;
// pub mod mmscanner;
pub mod omp_hack;
pub mod seqreader;
pub mod taxonomy;
pub mod utilities;

//  lookup accession numbers
// pub mod lookup_accession_numbers;
pub mod readcounts;
// pub mod mmap_file;
// pub mod omp_hack;
// pub mod utilities;

// k2mask
// pub mod k2mask;
// pub mod seqreader;

fn main() {
    println!("Hello, world!");
    // let nodes_content = "1\t|\t1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n2\t|\t1\t|\t2\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n";
    // let names_content = "1\t|\tall\t|\t\t|\tsynonym\t|\n1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\tBacteria <prokaryote>\t|\tscientific name\t|\n2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n2\t|\tProcaryotae\t|\tProcaryotae <Bacteria>\t|\tin-part\t|\n";

    // let nodes_file = "nodes.dmp";
    // let names_file = "names.dmp";

    // let mut file1 = File::create(nodes_file).unwrap();
    // let mut file2 = File::create(names_file).unwrap();

    // file1.write_all(nodes_content.as_bytes()).unwrap();
    // file2.write_all(names_content.as_bytes()).unwrap();

    let nodes_file = "/Users/mladenrasic/taxdump/nodes.dmp";
    let names_file = "/Users/mladenrasic/taxdump/names.dmp";

    let taxonomy = taxonomy::NCBITaxonomy::new(nodes_file, names_file).unwrap();

    // Serialize NCBITaxonomy to Kraken taxonomy and write into kraken.dmp
    taxonomy.convert_to_kraken_taxonomy("kraken.dmp").unwrap(); // writes to kraken_file

    // Deserialize
    let kraken_taxonomy = taxonomy::Taxonomy::from_file("kraken.dmp").unwrap();

    println!("{:#?}", kraken_taxonomy);

    assert_eq!(kraken_taxonomy.node_count, 2);
    assert_eq!(kraken_taxonomy.name_data_len, 12);
    assert_eq!(kraken_taxonomy.rank_data_len, 19);

    std::fs::remove_file(nodes_file).unwrap();
    std::fs::remove_file(names_file).unwrap();
    std::fs::remove_file("kraken.dmp").unwrap();
}
