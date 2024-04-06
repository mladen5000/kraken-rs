/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */


use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader, Write};

/// Taxonomy node structure, as stored in the Kraken 2 database.

#[derive(Clone, Debug)]
/// TaxonomyNode is a structure that holds the taxonomy node information.
/// A taxonomy node is a node in the NCBI taxonomy tree, containing information about the node's parent, children, name, rank, external ID, and godparent.
/// A node is a part of the NCBI taxonomy tree, which is a hierarchical classification of organisms.
/// In other words, a node is a taxonomic unit in the NCBI taxonomy tree, such as a species, genus, family, order, class, phylum, or kingdom.
///

#[derive(PartialEq, Serialize, Deserialize)]
pub struct TaxonomyNode {
    // Parent node
    pub parent_id: u64,
    // First child node
    pub first_child: u64,
    // Number of children
    pub child_count: u64,
    // idk
    pub name_offset: u64,
    // idk
    pub rank_offset: u64,
    // idk
    pub external_id: u64,
    // idk
    pub godparent_id: u64,
}

/// NCBITaxonomy is a structure that holds the NCBI taxonomy tree in memory.
/// Each map in the NCBITaxonomy structure is a mapping from a node ID to some information about the node.
/// A node ID is a unique identifier for a node in the NCBI taxonomy tree, obtained from the NCBI taxonomy dump files.
#[derive(Debug)]
pub struct NCBITaxonomy {
    /// Maps a node to its parent node.
    parent_map: BTreeMap<u64, u64>,
    /// Maps a node to its' scientific name
    name_map: BTreeMap<u64, String>,
    // Maps a node to its' rank
    rank_map: BTreeMap<u64, String>,
    /// Child_map is a map of parent IDs to child IDs.
    child_map: BTreeMap<u64, BTreeSet<u64>>,
    /// Marked_nodes is a set of node IDs that have been marked.
    /// In other words, a node is marked as "marked" if it is a descendant of a marked node.
    /// This is so that we can easily find all the descendants of a marked node.
    /// For example, if we mark the node for "Bacteria", then all the descendants of "Bacteria" will be marked as well.
    /// We use this later to convert the NCBI taxonomy to the Kraken 2 taxonomy format.
    marked_nodes: BTreeSet<u64>,

    /// The set of ranks {no rank, superkingdom, kingdom, phylum, class, order, family, genus, species}.
    known_ranks: BTreeSet<String>,
}

/// Create a new NCBITaxonomy object from dmp files.
impl NCBITaxonomy {
    pub fn new(
        nodes_filename: &str,
        names_filename: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        // 0. Initialize Maps and Sets
        let mut parent_map = BTreeMap::new();
        let mut name_map: BTreeMap<u64, String> = BTreeMap::new();
        let mut rank_map = BTreeMap::new();
        let mut child_map: BTreeMap<u64, BTreeSet<u64>> = BTreeMap::new();
        let mut marked_nodes = BTreeSet::new();
        let mut known_ranks = BTreeSet::new();

        // 1. Open files
        let nodes_file = BufReader::new(
            File::open(nodes_filename).map_err(|_| format!("Error opening {}", nodes_filename))?,
        );
        let names_file = BufReader::new(
            File::open(names_filename).map_err(|_| format!("Error opening {}", names_filename))?,
        );

        // 2. Process nodes
        // nodes.dmp format: <node_id | parent_id | rank | other | stuff | goes | here >
        for nodes_line in nodes_file.lines() {
            // Preprocess line
            let nodes_line = nodes_line.unwrap();
            let mut fields = nodes_line.split("\t|\t");

            // Extract fields:
            let node_id = fields.next().unwrap().parse::<u64>().unwrap();
            let mut parent_id = fields.next().unwrap().parse::<u64>().unwrap();
            parent_id = match node_id {
                0 => return Err("Bad parent id".into()),
                1 => 0,
                _ => parent_id,
            };
            let rank = fields.next().unwrap().to_string();

            // Build maps and sets
            parent_map.insert(node_id, parent_id);
            child_map.entry(parent_id).or_default().insert(node_id); // {Parent ID, {Node IDs}}
            rank_map.insert(node_id, rank.clone()); // {Node ID: Rank}
            known_ranks.insert(rank); // {Rank}
        }

        // 3. Process names.dmp
        // names.dmp format: <node_id | name | unique_name | name_type>
        for names_line in names_file.lines() {
            // Preprocess line
            let names_line = names_line.unwrap();
            let mut fields = names_line.split("\t|\t");

            // Parse line into fields
            let node_id = fields.next().unwrap().parse::<u64>().unwrap();
            let name = fields.next().unwrap().to_string();
            let _unique_name = fields.next();
            let name_type = fields.next().unwrap().trim_end_matches("\t|").to_string();

            // Build name_map
            if name_type == "scientific name" {
                name_map.insert(node_id, name);
            }
        }

        // 4. Idk but it's important
        marked_nodes.insert(1);

        // 5. Return the NCBITaxonomy object
        Ok(NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            marked_nodes,
            known_ranks,
        })
    }

    /// Mark a node and all its ancestors as "marked".
    pub fn mark_node(&mut self, taxid: u64) {
        let mut current_taxid = taxid;
        // While set of marked nodes does not contain the current taxid...
        while !self.marked_nodes.contains(&current_taxid) {
            // Add the current_taxid into the marked_nodes
            self.marked_nodes.insert(current_taxid);
            // Grab the parent node of the current_taxid
            current_taxid = *self.parent_map.get(&current_taxid).unwrap();
        }
    }

    /// Convert the NCBI taxonomy to the Kraken 2 taxonomy format.
    pub fn convert_to_kraken_taxonomy(
        &self,
        filename: &str,
    ) -> Result<Taxonomy, Box<dyn std::error::Error>> {
        let zeroes_node = TaxonomyNode {
            parent_id: 0,
            first_child: 0,
            child_count: 0,
            name_offset: 0,
            rank_offset: 0,
            external_id: 0,
            godparent_id: 0,
        };

        let mut taxo = Taxonomy::new();
        taxo.node_count = self.marked_nodes.len() + 1;
        taxo.nodes = vec![zeroes_node.clone(); taxo.node_count as usize];

        let mut rank_offsets = HashMap::new();
        let mut rank_data = String::new();

        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.push_str(rank);
            rank_data.push('\0');
        }

        let mut internal_node_id = 1;
        let mut external_id_map = BTreeMap::new();
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1);

        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(1);

        let mut name_data = String::new();

        while !bfs_queue.is_empty() {
            let external_node_id = bfs_queue.pop_front().unwrap();
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = zeroes_node.clone();
            node.parent_id = *external_id_map
                .get(&self.parent_map[&external_node_id])
                .ok_or_else(|| format!("Parent not found for node: {}", external_node_id))?;
            node.external_id = external_node_id;
            node.rank_offset = *rank_offsets
                .get(&self.rank_map[&external_node_id])
                .ok_or_else(|| format!("Rank not found for node: {}", external_node_id))?;
            node.name_offset = name_data.len() as u64;
            node.first_child = internal_node_id + 1 + bfs_queue.len() as u64;

            for &child_node in self
                .child_map
                .get(&external_node_id)
                .unwrap_or(&BTreeSet::new())
            {
                if self.marked_nodes.contains(&child_node) {
                    bfs_queue.push_back(child_node);
                    node.child_count += 1;
                }
            }

            taxo.nodes[internal_node_id as usize] = node;
            internal_node_id += 1;

            let blank_name = String::new();
            let name = self.name_map.get(&external_node_id).unwrap_or(&blank_name);
            name_data.push_str(name);
            name_data.push('\0');
        }

        taxo.rank_data = rank_data.into_bytes();
        taxo.rank_data_len = taxo.rank_data.len();

        taxo.name_data = name_data.into_bytes();
        taxo.name_data_len = taxo.name_data.len();

        taxo.write_to_disk(filename)?;

        Ok(taxo)
    }
}

/// Taxonomy is a structure that holds the taxonomy tree in memory.
/// It differs from NCBITaxonomy and TaxonomyNode by storing the taxonomy tree in a file-backed format.
#[derive(Debug, Serialize, Deserialize)]
pub struct Taxonomy {
    /// File_backed is a boolean that indicates whether the taxonomy is file-backed.
    pub file_backed: bool,
    /// A vector of taxonomy nodes.
    pub nodes: Vec<TaxonomyNode>,
    /// The number of nodes in the taxonomy.
    pub node_count: usize,
    /// A vector of name data.
    pub name_data: Vec<u8>,
    /// The length of the name data.
    pub name_data_len: usize,
    /// A vector of rank data.
    pub rank_data: Vec<u8>,
    /// The length of the rank data.
    pub rank_data_len: usize,
    /// A map of external IDs to internal IDs.
    pub external_to_internal_id_map: HashMap<u64, u64>,
}

impl Taxonomy {
    pub fn new() -> Self {
        Taxonomy {
            file_backed: false,
            nodes: Vec::new(),
            node_count: 0,
            name_data: Vec::new(),
            name_data_len: 0,
            rank_data: Vec::new(),
            rank_data_len: 0,
            external_to_internal_id_map: HashMap::new(),
        }
    }

    // pub fn from_file_old(filename: &str, memory_mapping: bool) -> Self {
    //     let mut taxo = Taxonomy::new();
    //     taxo.init(filename, memory_mapping);
    //     taxo
    // }

    /// Majority of this code is building the Taxonomy struct from a file.
    /// Each struct field is a buffer which we fill from the disk file.
    /// Basically a deserialization function.
    fn init(&mut self, filename: &str, memory_mapping: bool) {
        if memory_mapping {
            // TODO: Implement memory mapping
        } else {
            let mut file = BufReader::new(File::open(filename).unwrap());
            let _magic = [0; 8]; // 8 bytes possibly k2taxdat?

            // ...

            // std::slice::from_mut(&mut magic) takes in a &mut u8 and returns a &mut [u8]
            // Read exactly 8 bytes from the file into the magic buffer
            file.read_exact(std::slice::from_mut(&mut (self.node_count as u8)))
                .unwrap();
            file.read_exact(std::slice::from_mut(&mut (self.name_data_len as u8)))
                .unwrap();
            file.read_exact(std::slice::from_mut(&mut (self.rank_data_len as u8)))
                .unwrap();

            // Initialize a Vec<taxnodes> with size node_count
            self.nodes = vec![
                TaxonomyNode {
                    parent_id: 0,
                    first_child: 0,
                    child_count: 0,
                    name_offset: 0,
                    rank_offset: 0,
                    external_id: 0,
                    godparent_id: 0,
                };
                self.node_count
            ];

            // Read exactly (node_count * size_of::<TaxonomyNode>) bytes from the file into the nodes vector
            // Go into file, read exactly (node_count * size_of::<TaxonomyNode>) bytes, and put it into the nodes vector
            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    // Build a slice from ptr + length
                    self.nodes.as_mut_ptr() as *mut u8, // Return a mut pointer to the first element of the slice
                    self.node_count * std::mem::size_of::<TaxonomyNode>(), // Return the size of the slice in bytes
                )
            })
            .unwrap();

            self.name_data = vec![0; self.name_data_len];
            file.read_exact(&mut self.name_data).unwrap();

            self.rank_data = vec![0; self.rank_data_len];
            file.read_exact(&mut self.rank_data).unwrap();
        }
    }

    /// Check if a node is an ancestor of another node.
    pub fn is_a_ancestor_of_b(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        let mut current_b = b;
        while current_b > a {
            current_b = self.nodes[current_b as usize].parent_id;
        }
        current_b == a
    }

    /// Find the lowest common ancestor of two nodes.
    pub fn lowest_common_ancestor(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }
        let mut current_a = a;
        let mut current_b = b;
        // If the nodes differ
        println!("Current A: {} Current B: {}", current_a, current_b);
        while current_a != current_b {
            // Get the parent of the larger one and try again with its parent
            if current_a > current_b {
                current_a = self.nodes[current_a as usize].parent_id;
            } else {
                current_b = self.nodes[current_b as usize].parent_id;
            }
        }
        current_a
    }

    /// Serialize the taxonomy to disk.
    // pub fn write_to_disk_old(&self, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    //     let mut file = File::create(filename)?;
    //     file.write_all(b"K2TAXDAT")?;
    //     file.write_all(&self.node_count.to_le_bytes())?;
    //     file.write_all(&self.name_data_len.to_le_bytes())?;
    //     file.write_all(&self.rank_data_len.to_le_bytes())?;
    //     file.write_all(unsafe {
    //         std::slice::from_raw_parts(
    //             self.nodes.as_ptr() as *const u8,
    //             self.node_count * std::mem::size_of::<TaxonomyNode>(),
    //         )
    //     })?;
    //     file.write_all(&self.name_data)?;
    //     file.write_all(&self.rank_data)?;
    //     Ok(())
    // }

    /// Generate a mapping from external IDs to internal IDs.   
    pub fn generate_external_to_internal_id_map(&mut self) -> &mut Self {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);
        for i in 1..self.node_count {
            self.external_to_internal_id_map
                .insert(self.nodes[i].external_id, i as u64);
        }
        self
    }

    /// Get the internal ID of a node given its external ID.
    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        *self
            .external_to_internal_id_map
            .get(&external_id)
            .unwrap_or(&0)
    }
    pub fn write_to_disk(&self, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(filename)?;
        let encoded: Vec<u8> = serde_json::to_vec(self)?;
        file.write_all(&encoded)?;
        Ok(())
    }

    pub fn from_file(filename: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = File::open(filename)?;
        let mut encoded = Vec::new();
        file.read_to_end(&mut encoded)?;
        let taxo: Taxonomy = serde_json::from_slice(&encoded)?;
        Ok(taxo)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_ncbi_taxonomy_creation() {
        // Ok
        // Fake data..
        let nodes_content = "1\t|\t1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n2\t|\t1\t|\t2\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n";
        let names_content = "1\t|\tall\t|\t\t|\tsynonym\t|\n1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\tBacteria <prokaryote>\t|\tscientific name\t|\n2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n2\t|\tProcaryotae\t|\tProcaryotae <Bacteria>\t|\tin-part\t|\n";

        let nodes_file = "nodes.dmp";
        let names_file = "names.dmp";

        let mut file1 = File::create(nodes_file).unwrap();
        let mut file2 = File::create(names_file).unwrap();

        file1.write_all(nodes_content.as_bytes()).unwrap();
        file2.write_all(names_content.as_bytes()).unwrap();

        // Or...
        // let nodes_file = "/Users/mladenrasic/taxdump/nodes.dmp";
        // let names_file = "/Users/mladenrasic/taxdump/names.dmp";

        let ncbitaxonomy = NCBITaxonomy::new(nodes_file, names_file).unwrap();
        // println!("{:#?}", taxonomy);

        assert_eq!(ncbitaxonomy.parent_map.len(), 2);
        assert_eq!(ncbitaxonomy.name_map.len(), 2);
        assert_eq!(ncbitaxonomy.rank_map.len(), 2);
        assert_eq!(ncbitaxonomy.child_map.len(), 2);
        assert_eq!(ncbitaxonomy.marked_nodes.len(), 1);
        assert_eq!(ncbitaxonomy.known_ranks.len(), 2);

        std::fs::remove_file(nodes_file).unwrap();
        std::fs::remove_file(names_file).unwrap();
    }

    #[test]
    fn test_ncbi_taxonomy_mark_node() {
        // Ok
        let nodes_content = "1\t|\t1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n2\t|\t1\t|\t2\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n";
        let names_content = "1\t|\tall\t|\t\t|\tsynonym\t|\n1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\tBacteria <prokaryote>\t|\tscientific name\t|\n2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n2\t|\tProcaryotae\t|\tProcaryotae <Bacteria>\t|\tin-part\t|\n";

        let nodes_file = "nodes.dmp";
        let names_file = "names.dmp";

        let mut file1 = File::create(nodes_file).unwrap();
        let mut file2 = File::create(names_file).unwrap();

        file1.write_all(nodes_content.as_bytes()).unwrap();
        file2.write_all(names_content.as_bytes()).unwrap();

        let mut ncbitaxonomy = NCBITaxonomy::new(nodes_file, names_file).unwrap();

        ncbitaxonomy.mark_node(2);

        assert_eq!(ncbitaxonomy.marked_nodes.len(), 2);
        assert!(ncbitaxonomy.marked_nodes.contains(&1));
        assert!(ncbitaxonomy.marked_nodes.contains(&2));

        std::fs::remove_file(nodes_file).unwrap();
        std::fs::remove_file(names_file).unwrap();
    }

    #[test]
    fn test_ncbi_taxonomy_convert_to_kraken_taxonomy() {
        let nodes_content = "1\t|\t1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n2\t|\t1\t|\t2\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n";
        let names_content = "1\t|\tall\t|\t\t|\tsynonym\t|\n1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\tBacteria <prokaryote>\t|\tscientific name\t|\n2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n2\t|\tProcaryotae\t|\tProcaryotae <Bacteria>\t|\tin-part\t|\n";

        let nodes_file = "nodes.dmp";
        let names_file = "names.dmp";

        let mut file1 = File::create(nodes_file).unwrap();
        let mut file2 = File::create(names_file).unwrap();

        file1.write_all(nodes_content.as_bytes()).unwrap();
        file2.write_all(names_content.as_bytes()).unwrap();

        // let nodes_file = "/Users/mladenrasic/taxdump/nodes.dmp";
        // let names_file = "/Users/mladenrasic/taxdump/names.dmp";

        let taxonomy = NCBITaxonomy::new(nodes_file, names_file).unwrap();

        // Serialize NCBITaxonomy to Kraken taxonomy and write into kraken.dmp
        taxonomy.convert_to_kraken_taxonomy("kraken.dmp").unwrap(); // writes to kraken_file

        // Deserialize //error here
        let kraken_taxonomy = Taxonomy::from_file("kraken.dmp").unwrap();

        println!("{:#?}", kraken_taxonomy);

        assert_eq!(kraken_taxonomy.node_count, 2);
        assert_eq!(kraken_taxonomy.name_data_len, 5);
        assert_eq!(kraken_taxonomy.rank_data_len, 4);

        // std::fs::remove_file(nodes_file).unwrap();
        // std::fs::remove_file(names_file).unwrap();
        // std::fs::remove_file("kraken.dmp").unwrap();
    }

    #[test]
    fn test_taxonomy_is_a_ancestor_of_b() {
        let nodes_content = "1\t|\t1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n2\t|\t1\t|\t2\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n3\t|\t1\t|\t3\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n";
        let names_content = "1\t|\tall\t|\t\t|\tsynonym\t|\n1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\tBacteria <prokaryote>\t|\tscientific name\t|\n2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n2\t|\tProcaryotae\t|\tProcaryotae <Bacteria>\t|\tin-part\t|\n3\t|\tFirmicutes\t|\tFirmicutes <Bacteria>\t|\tscientific name\t|\n3\t|\tActinobacteria\t|\tActinobacteria <Bacteria>\t|\tscientific name\t|\n";

        let nodes_file = "nodes.dmp";
        let names_file = "names.dmp";

        let mut file = File::create(nodes_file).unwrap();
        file.write_all(nodes_content.as_bytes()).unwrap();

        let mut file = File::create(names_file).unwrap();
        file.write_all(names_content.as_bytes()).unwrap();

        // let nodes_file = "/Users/mladenrasic/taxdump/nodes.dmp";
        // let names_file = "/Users/mladenrasic/taxdump/names.dmp";

        let ncbi_taxonomy = NCBITaxonomy::new(nodes_file, names_file).unwrap();
        let mut kraken_taxonomy = ncbi_taxonomy
            .convert_to_kraken_taxonomy("kraken.dmp")
            .unwrap();

        let mut_taxonomy = kraken_taxonomy.generate_external_to_internal_id_map();

        println!("{:#?}", mut_taxonomy);

        assert_eq!(mut_taxonomy.external_to_internal_id_map.len(), 2);
        assert_eq!(mut_taxonomy.external_to_internal_id_map.get(&1), Some(&1));
        assert_eq!(mut_taxonomy.external_to_internal_id_map.get(&0), Some(&0));
    }
}
