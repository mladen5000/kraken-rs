use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, Read};
use std::io::{BufWriter, Write};
use std::path::Path;
const FILE_MAGIC: &str = "K2TAXDAT";

#[derive(Debug)]
pub struct NCBITaxonomy {
    pub parent_map: HashMap<u64, u64>,
    pub name_map: HashMap<u64, String>,
    pub rank_map: HashMap<u64, String>,
    pub child_map: HashMap<u64, HashSet<u64>>,
    pub marked_nodes: HashSet<u64>,
    pub known_ranks: HashSet<String>,
}

impl NCBITaxonomy {
    /// Creates a new NCBITaxonomy by reading nodes and names files.
    ///
    /// # Arguments
    ///
    /// * `nodes_filename` - The path to the nodes file.
    /// * `names_filename` - The path to the names file.
    pub fn new(nodes_filename: &str, names_filename: &str) -> io::Result<Self> {
        let mut taxonomy = NCBITaxonomy {
            parent_map: HashMap::new(),
            name_map: HashMap::new(),
            rank_map: HashMap::new(),
            child_map: HashMap::new(),
            marked_nodes: HashSet::new(),
            known_ranks: HashSet::new(),
        };

        taxonomy.read_nodes(nodes_filename)?;
        taxonomy.read_names(names_filename)?;
        Ok(taxonomy)
    }

    fn read_nodes(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?.trim_end_matches("\t|").to_string();
            let tokens: Vec<&str> = line.split("\t|\t").collect();
            let node_id = tokens[0].parse::<u64>().unwrap();
            let parent_id = tokens[1].parse::<u64>().unwrap();
            let rank = tokens[2].to_string();

            if node_id == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "node ID == 0"));
            }

            self.parent_map.insert(node_id, parent_id);
            self.child_map
                .entry(parent_id)
                .or_insert_with(HashSet::new)
                .insert(node_id);
            self.rank_map.insert(node_id, rank.clone());
            self.known_ranks.insert(rank.clone());
        }

        Ok(())
    }

    fn read_names(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?.trim_end_matches("\t|").to_string();
            let tokens: Vec<&str> = line.split("\t|\t").collect();
            let node_id = tokens[0].parse::<u64>().unwrap();
            let name = tokens[1].to_string();

            if tokens[3] == "scientific name" {
                self.name_map.insert(node_id, name);
            }
        }

        Ok(())
    }

    pub fn mark_node(&mut self, mut taxid: u64) {
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            taxid = *self.parent_map.get(&taxid).unwrap_or(&0);
        }
    }

    /// Converts the NCBI taxonomy to the Kraken2 taxonomy format and writes it to disk.
    ///
    /// # Arguments
    ///
    /// * `filename` - The path to the output file.
    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> io::Result<()> {
        let mut taxo_nodes: Vec<TaxonomyNode> = Vec::new();
        let mut name_data = String::new();
        let mut rank_data = String::new();
        let mut rank_offsets = HashMap::new();

        // Collect all unique ranks and prepare rank data with offsets
        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.push_str(rank);
            rank_data.push('\0'); // Null-terminated strings
        }

        let mut bfs_queue = VecDeque::new();
        let mut external_to_internal_id_map = HashMap::new();

        bfs_queue.push_back(1); // Start BFS from the root node, which is usually '1'
        external_to_internal_id_map.insert(1, 1);

        while let Some(external_id) = bfs_queue.pop_front() {
            let internal_id = external_to_internal_id_map[&external_id];
            let parent_id = self.parent_map[&external_id];
            let internal_parent_id = *external_to_internal_id_map.get(&parent_id).unwrap_or(&0);

            let node = TaxonomyNode {
                parent_id: internal_parent_id,
                external_id,
                first_child: internal_id + bfs_queue.len() as u64 + 1, // Estimate first_child index
                child_count: 0,                                        // Will count children below
                name_offset: name_data.len() as u64,
                rank_offset: *rank_offsets.get(&self.rank_map[&external_id]).unwrap(),
                godparent_id: 0, // Not used in this simplified version
            };

            // Add node name to the name data pool
            if let Some(name) = self.name_map.get(&external_id) {
                name_data.push_str(name);
                name_data.push('\0'); // Null-terminated strings
            }

            // Queue children for BFS
            if let Some(children) = self.child_map.get(&external_id) {
                for &child_id in children {
                    if self.marked_nodes.contains(&child_id) {
                        bfs_queue.push_back(child_id);
                        external_to_internal_id_map
                            .insert(child_id, internal_id + bfs_queue.len() as u64);
                        taxo_nodes[internal_id as usize - 1].child_count += 1;
                    }
                }
            }

            taxo_nodes.push(node);
        }

        // Writing to disk
        let mut file = BufWriter::new(File::create(Path::new(filename))?);
        file.write_all(b"K2TAXDAT")?;
        for node in &taxo_nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            file.write_all(&node.first_child.to_le_bytes())?;
            file.write_all(&node.child_count.to_le_bytes())?;
            file.write_all(&node.name_offset.to_le_bytes())?;
            file.write_all(&node.rank_offset.to_le_bytes())?;
            file.write_all(&node.external_id.to_le_bytes())?;
            file.write_all(&node.godparent_id.to_le_bytes())?;
        }
        file.write_all(name_data.as_bytes())?;
        file.write_all(rank_data.as_bytes())?;

        Ok(())
    }
}
#[derive(Debug, Clone)]
pub struct Taxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: Vec<u8>,
    pub rank_data: Vec<u8>,
    node_count: usize,
    file_backed: bool,
}

impl Taxonomy {
    pub fn new(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        let mut taxonomy = Taxonomy {
            nodes: Vec::new(),
            name_data: Vec::new(),
            rank_data: Vec::new(),
            node_count: 0,
            file_backed: memory_mapping,
        };

        if memory_mapping {
            // Implement memory mapping logic here. Rust doesn't directly support it in the standard library.
            // External crates like `memmap` or `mmap` might be necessary.
            // For demonstration, this will be left as a conceptual placeholder.
            println!("Memory mapping is currently not implemented.");
        } else {
            taxonomy.init_from_file(filename)?;
        }

        Ok(taxonomy)
    }

    fn init_from_file(&mut self, filename: &str) -> io::Result<()> {
        let file = fs::read(filename)?;
        let mut offset = 0;

        // Check for file magic
        let magic = &file[offset..offset + FILE_MAGIC.len()];
        offset += FILE_MAGIC.len();
        if magic != FILE_MAGIC.as_bytes() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed taxonomy file",
            ));
        }

        // Read node count
        self.node_count = u64::from_le_bytes(file[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;

        // Read name data length
        let name_data_len =
            u64::from_le_bytes(file[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;

        // Read rank data length
        let rank_data_len =
            u64::from_le_bytes(file[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;

        // Read nodes
        self.nodes = Vec::with_capacity(self.node_count);
        for _ in 0..self.node_count {
            let node = TaxonomyNode {
                parent_id: u64::from_le_bytes(file[offset..offset + 8].try_into().unwrap()),
                first_child: u64::from_le_bytes(file[offset + 8..offset + 16].try_into().unwrap()),
                child_count: u64::from_le_bytes(file[offset + 16..offset + 24].try_into().unwrap()),
                name_offset: u64::from_le_bytes(file[offset + 24..offset + 32].try_into().unwrap()),
                rank_offset: u64::from_le_bytes(file[offset + 32..offset + 40].try_into().unwrap()),
                external_id: u64::from_le_bytes(file[offset + 40..offset + 48].try_into().unwrap()),
                godparent_id: u64::from_le_bytes(
                    file[offset + 48..offset + 56].try_into().unwrap(),
                ),
            };
            self.nodes.push(node);
            offset += 56; // size of one TaxonomyNode
        }

        // Read name and rank data
        self.name_data = file[offset..offset + name_data_len].to_vec();
        offset += name_data_len;
        self.rank_data = file[offset..offset + rank_data_len].to_vec();

        Ok(())
    }

    pub fn write_to_disk(&self, filename: &str) -> io::Result<()> {
        let mut file = BufWriter::new(OpenOptions::new().write(true).create(true).open(filename)?);

        file.write_all(FILE_MAGIC.as_bytes())?;
        file.write_all(&(self.node_count as u64).to_le_bytes())?;
        file.write_all(&(self.name_data.len() as u64).to_le_bytes())?;
        file.write_all(&(self.rank_data.len() as u64).to_le_bytes())?;
        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;
        for node in &self.nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            file.write_all(&node.first_child.to_le_bytes())?;
            file.write_all(&node.child_count.to_le_bytes())?;
            file.write_all(&node.name_offset.to_le_bytes())?;
            file.write_all(&node.rank_offset.to_le_bytes())?;
            file.write_all(&node.external_id.to_le_bytes())?;
            file.write_all(&node.godparent_id.to_le_bytes())?;
        }

        Ok(())
    }
}

use std::convert::TryInto;

#[derive(Clone, Debug)]
pub struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64, // Reserved for potential future optimizations
}

impl TaxonomyNode {
    /// Creates a `TaxonomyNode` from a slice of bytes.
    /// This function assumes the bytes are in the correct order and format.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let parent_id = u64::from_le_bytes(bytes[0..8].try_into().unwrap());
        let first_child = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        let child_count = u64::from_le_bytes(bytes[16..24].try_into().unwrap());
        let name_offset = u64::from_le_bytes(bytes[24..32].try_into().unwrap());
        let rank_offset = u64::from_le_bytes(bytes[32..40].try_into().unwrap());
        let external_id = u64::from_le_bytes(bytes[40..48].try_into().unwrap());
        let godparent_id = u64::from_le_bytes(bytes[48..56].try_into().unwrap());

        TaxonomyNode {
            parent_id,
            first_child,
            child_count,
            name_offset,
            rank_offset,
            external_id,
            godparent_id,
        }
    }
}

// Add any additional methods such as IsAAncestorOfB and LowestCommonAncestor here.

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::io::Cursor;

//     const NODES_DATA: &str =
//         "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
// 2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
// 3\t|\t1\t|\tclade\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|";

//     const NAMES_DATA: &str = "1\t|\tall\t|\t\t|\tsynonym\t|
// 1\t|\troot\t|\t\t|\tscientific name\t|
// 2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|
// 3\t|\tAT-rich\t|\t\t|\tsynonym\t|";

//     #[test]
//     fn test_ncbi_taxonomy_new() {
//         let nodes_file = Cursor::new(NODES_DATA);
//         let names_file = Cursor::new(NAMES_DATA);

//         let taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         println!("{:#?}", taxonomy);

//         assert_eq!(taxonomy.parent_map.len(), 3);
//         assert_eq!(taxonomy.name_map.len(), 2);
//         assert_eq!(taxonomy.rank_map.len(), 3);
//         assert_eq!(taxonomy.child_map.len(), 1);
//         assert_eq!(taxonomy.marked_nodes.len(), 1);
//         assert_eq!(taxonomy.known_ranks.len(), 3);

//         assert_eq!(taxonomy.parent_map[&1], 0);
//         assert_eq!(taxonomy.parent_map[&2], 1);
//         assert_eq!(taxonomy.parent_map[&3], 1);

//         assert_eq!(taxonomy.name_map[&1], "root");
//         assert_eq!(taxonomy.name_map[&2], "Bacteria");
//         // assert_eq!(taxonomy.name_map[&3], "AT-rich"); // Not "scientific name"

//         assert_eq!(taxonomy.rank_map[&1], "no rank");
//         assert_eq!(taxonomy.rank_map[&2], "superkingdom");
//         assert_eq!(taxonomy.rank_map[&3], "clade");

//         assert!(taxonomy.child_map[&1].contains(&2));
//         assert!(taxonomy.child_map[&1].contains(&3));

//         assert!(taxonomy.marked_nodes.contains(&1));

//         assert!(taxonomy.known_ranks.contains("no rank"));
//         assert!(taxonomy.known_ranks.contains("superkingdom"));
//         assert!(taxonomy.known_ranks.contains("clade"));
//     }

//     #[test]
//     fn test_ncbi_taxonomy_mark_node() {
//         let nodes_file = Cursor::new(NODES_DATA);
//         let names_file = Cursor::new(NAMES_DATA);

//         let mut taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         println!("{:?}", taxonomy);

//         taxonomy.mark_node(3);

//         assert!(taxonomy.marked_nodes.contains(&1));
//         assert!(taxonomy.marked_nodes.contains(&3));
//     }

//     #[test]
//     fn test_taxonomy_is_a_ancestor_of_b() {
//         let nodes_file = Cursor::new(NODES_DATA);
//         let names_file = Cursor::new(NAMES_DATA);

//         let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

//         let mut taxonomy = Taxonomy::new("test_taxonomy.k2d");
//         taxonomy.generate_external_to_internal_id_map();

//         let root_id = taxonomy.get_internal_id(1);
//         let bacteria_id = taxonomy.get_internal_id(2);
//         let at_rich_id = taxonomy.get_internal_id(3);

//         assert!(taxonomy.is_a_ancestor_of_b(root_id, bacteria_id));
//         assert!(taxonomy.is_a_ancestor_of_b(root_id, at_rich_id));
//         assert!(!taxonomy.is_a_ancestor_of_b(bacteria_id, root_id));
//         assert!(!taxonomy.is_a_ancestor_of_b(bacteria_id, at_rich_id));
//         assert!(!taxonomy.is_a_ancestor_of_b(at_rich_id, root_id));
//         assert!(!taxonomy.is_a_ancestor_of_b(at_rich_id, bacteria_id));
//     }

//     #[test]
//     fn test_taxonomy_lowest_common_ancestor() {
//         let nodes_file = Cursor::new(NODES_DATA);
//         let names_file = Cursor::new(NAMES_DATA);

//         let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

//         let taxonomy = Taxonomy::new("test_taxonomy.k2d");

//         assert_eq!(taxonomy.lowest_common_ancestor(1, 2), 1);
//         assert_eq!(taxonomy.lowest_common_ancestor(1, 3), 1);
//         assert_eq!(taxonomy.lowest_common_ancestor(2, 3), 1);
//         assert_eq!(taxonomy.lowest_common_ancestor(2, 2), 2);
//         assert_eq!(taxonomy.lowest_common_ancestor(3, 3), 3);
//     }

//     #[test]
//     fn test_taxonomy_get_internal_id() {
//         let nodes_file = Cursor::new(NODES_DATA);
//         let names_file = Cursor::new(NAMES_DATA);

//         let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

//         let mut taxonomy = Taxonomy::new("test_taxonomy.k2d");
//         taxonomy.generate_external_to_internal_id_map();

//         assert_eq!(taxonomy.get_internal_id(1), 1);
//         assert_eq!(taxonomy.get_internal_id(2), 2);
//         assert_eq!(taxonomy.get_internal_id(3), 3);
//         assert_eq!(taxonomy.get_internal_id(4), 0);
//     }

//     #[cfg(test)]
//     // Existing constants for nodes and names data.

//     // Test handling of empty input files.
//     #[test]
//     fn test_empty_input_files() {
//         let nodes_file = Cursor::new("");
//         let names_file = Cursor::new("");

//         let taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
//         assert!(taxonomy.parent_map.is_empty());
//         assert!(taxonomy.name_map.is_empty());
//         assert!(taxonomy.rank_map.is_empty());
//         assert!(taxonomy.child_map.is_empty());
//         assert!(taxonomy.marked_nodes.is_empty());
//         assert!(taxonomy.known_ranks.is_empty());
//     }

//     // Test handling of input files with malformed lines.
//     #[test]
//     fn test_malformed_lines() {
//         let nodes_data =
//             "1\t|\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|";
//         let names_data = "1\t|\troot\t|\t\t|\t";
//         let nodes_file = Cursor::new(nodes_data);
//         let names_file = Cursor::new(names_data);

//         let result = std::panic::catch_unwind(|| {
//             NCBITaxonomy::from_readers(nodes_file, names_file);
//         });

//         assert!(result.is_err(), "Expected a panic due to malformed input");
//     }

//     // Test handling of lines with missing fields.
//     #[test]
//     fn test_incomplete_lines() {
//         let nodes_data = "1\t|\t1\t|\tno rank";
//         let names_data = "1\t|\troot";
//         let nodes_file = Cursor::new(nodes_data);
//         let names_file = Cursor::new(names_data);

//         let result = std::panic::catch_unwind(|| {
//             NCBITaxonomy::from_readers(nodes_file, names_file);
//         });

//         assert!(result.is_err(), "Expected a panic due to incomplete input");
//     }

//     // Test correct handling and parsing of node and name types.
//     #[test]
//     fn test_valid_input_handling() {
//         let nodes_data =
//             "2\t|\t1\t|\tgenus\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|";
//         let names_data = "2\t|\tEscherichia coli\t|\t\t|\tscientific name\t|";
//         let nodes_file = Cursor::new(nodes_data);
//         let names_file = Cursor::new(names_data);

//         let taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);

//         assert_eq!(taxonomy.parent_map.get(&2), Some(&1));
//         assert_eq!(
//             taxonomy.name_map.get(&2),
//             Some(&"Escherichia coli".to_string())
//         );
//         assert_eq!(taxonomy.rank_map.get(&2), Some(&"genus".to_string()));
//         assert!(taxonomy.child_map[&1].contains(&2));
//     }

//     // Test that nodes with non-scientific names are not included in the name_map.
//     #[test]
//     fn test_non_scientific_names_exclusion() {
//         let names_data = "3\t|\tE. coli\t|\t\t|\tcommon name\t|";
//         let names_file = Cursor::new(names_data);

//         let nodes_data =
//             "3\t|\t2\t|\tspecies\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|";
//         let nodes_file = Cursor::new(nodes_data);

//         let taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);

//         assert!(
//             taxonomy.name_map.get(&3).is_none(),
//             "Common names should not be included in the name_map."
//         );
//     }
// }
