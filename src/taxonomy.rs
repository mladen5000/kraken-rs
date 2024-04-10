use libc::if_nameindex;
use memmap::Mmap;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

const FILE_MAGIC: &[u8] = b"K2TAXDAT";

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
    pub fn new<P: AsRef<Path>>(nodes_filename: P, names_filename: P) -> Self {
        let nodes_file =
            BufReader::new(File::open(nodes_filename).expect("Failed to open nodes file"));
        let names_file =
            BufReader::new(File::open(names_filename).expect("Failed to open names file"));
        Self::from_readers(nodes_file, names_file)
    }
    fn from_readers<R: Read>(nodes_reader: R, names_reader: R) -> Self {
        let mut parent_map = HashMap::new();
        let mut name_map = HashMap::new();
        let mut rank_map = HashMap::new();
        let mut child_map = HashMap::new();
        let mut marked_nodes = HashSet::new();
        let mut known_ranks = HashSet::new();

        for line in BufReader::new(nodes_reader).lines() {
            let line = line.expect("Failed to read line from nodes file");
            let fields: Vec<&str> = line.trim().split("\t|\t").collect();

            let node_id = fields[0].parse().expect("Failed to parse node ID");
            let parent_id = fields[1].parse().unwrap_or(0);
            let rank = fields[2].to_string();

            if node_id == 0 {
                panic!("Attempt to create taxonomy with node ID 0");
            }

            if node_id == 1 {
                parent_map.insert(node_id, 0);
            } else {
                parent_map.insert(node_id, parent_id);
            }

            child_map
                .entry(parent_id)
                .or_insert_with(HashSet::new)
                .insert(node_id);
            rank_map.insert(node_id, rank.clone());
            known_ranks.insert(rank);
        }

        for line in BufReader::new(names_reader).lines() {
            let line = line.expect("Failed to read line from names file");
            let fields: Vec<&str> = line.trim().split("\t|\t").collect();

            let node_id = fields[0].parse().expect("Failed to parse node ID");
            let name = fields[1].to_string();
            let name_type = fields[3];

            let name_type = name_type
                .strip_suffix("\t|")
                .unwrap_or(name_type)
                .to_string();

            if node_id == 0 {
                panic!("Attempt to create taxonomy with node ID 0");
            }

            if name_type == "scientific name" {
                name_map.insert(node_id, name);
            }
        }

        marked_nodes.insert(1); // Mark root node

        NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            marked_nodes,
            known_ranks,
        }
    }

    /// Marks a given taxonomy node and all its unmarked ancestors.
    ///
    /// # Arguments
    ///
    /// * `taxid` - The taxonomy ID of the node to mark.
    pub fn mark_node(&mut self, taxid: u64) {
        let mut current_id = taxid;
        while !self.marked_nodes.contains(&current_id) {
            self.marked_nodes.insert(current_id);
            current_id = *self
                .parent_map
                .get(&current_id)
                .expect("Parent ID not found");
        }
    }

    /// Converts the NCBI taxonomy to the Kraken2 taxonomy format and writes it to disk.
    ///
    /// # Arguments
    ///
    /// * `filename` - The path to the output file.
    fn convert_to_kraken_taxonomy<P: AsRef<Path>>(&self, filename: P) {
        let mut nodes = Vec::new();
        let mut name_data = String::new();
        let mut rank_data = String::new();
        let mut rank_offsets = HashMap::new();
        let mut external_id_map = HashMap::new();

        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.push_str(rank);
            rank_data.push('\0');
        }

        let mut internal_node_id = 0;
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1);

        let mut bfs_queue = std::collections::VecDeque::new();
        bfs_queue.push_back(1);

        while let Some(external_node_id) = bfs_queue.pop_front() {
            internal_node_id += 1;
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = TaxonomyNode {
                parent_id: external_id_map[&self.parent_map[&external_node_id]],
                first_child: 0,
                child_count: 0,
                name_offset: name_data.len() as u64,
                rank_offset: rank_offsets[&self.rank_map[&external_node_id]],
                external_id: external_node_id,
                godparent_id: 0,
            };

            let name = match self.name_map.get(&external_node_id) {
                Some(name) => name.to_string(),
                None => format!("Unnamed taxon {}", external_node_id),
            };
            name_data.push_str(&name);
            name_data.push('\0');

            let mut child_nodes = Vec::new();
            for child_node in self
                .child_map
                .get(&external_node_id)
                .expect("Child nodes not found")
            {
                if self.marked_nodes.contains(child_node) {
                    child_nodes.push(*child_node);
                    node.child_count += 1;
                }
            }

            if !child_nodes.is_empty() {
                node.first_child = internal_node_id + 1;
            }

            nodes.push(node);
            bfs_queue.extend(child_nodes);
        }

        let mut taxonomy = Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map: HashMap::new(),
        };

        taxonomy.write_to_disk(filename);
    }
}
#[derive(Debug, Serialize, Deserialize)]
struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

#[derive(Debug)]
pub struct Taxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: String,
    pub rank_data: String,
    pub external_to_internal_id_map: HashMap<u64, u64>,
}

impl Taxonomy {
    /// Creates a new Taxonomy by reading from a file.
    ///
    /// # Arguments
    ///
    /// * `filename` - The path to the taxonomy file.
    pub fn new<P: AsRef<Path>>(filename: P) -> Self {
        let file = File::open(filename).expect("Failed to open taxonomy file");
        let mmap = unsafe { Mmap::map(&file).expect("Failed to memory-map taxonomy file") };

        let mut cursor = 0;
        let magic = &mmap[cursor..cursor + FILE_MAGIC.len()];
        assert_eq!(magic, FILE_MAGIC, "Invalid taxonomy file format");
        cursor += FILE_MAGIC.len();

        let node_count = Self::read_u64(&mmap, &mut cursor);
        let name_data_len = Self::read_u64(&mmap, &mut cursor);
        let rank_data_len = Self::read_u64(&mmap, &mut cursor);

        let nodes_size = node_count as usize * std::mem::size_of::<TaxonomyNode>();
        let nodes_data = &mmap[cursor..cursor + nodes_size];
        cursor += nodes_size;

        let nodes = bincode::deserialize(nodes_data).expect("Failed to deserialize taxonomy nodes");

        let name_data =
            String::from_utf8_lossy(&mmap[cursor..cursor + name_data_len as usize]).to_string();
        cursor += name_data_len as usize;

        let rank_data =
            String::from_utf8_lossy(&mmap[cursor..cursor + rank_data_len as usize]).to_string();

        Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map: HashMap::new(),
        }
    }

    /// Reads a u64 value from the memory-mapped file.
    ///
    /// # Arguments
    ///
    /// * `mmap` - The memory-mapped file.
    /// * `cursor` - The current cursor position.
    ///
    /// # Returns
    ///
    /// The read u64 value.
    pub fn read_u64(mmap: &[u8], cursor: &mut usize) -> u64 {
        let value =
            bincode::deserialize(&mmap[*cursor..*cursor + 8]).expect("Failed to read u64 value");
        *cursor += 8;
        value
    }

    /// Checks if node A is an ancestor of node B.
    ///
    /// # Arguments
    ///
    /// * `a` - The taxonomy ID of node A.
    /// * `b` - The taxonomy ID of node B.
    ///
    /// # Returns
    ///
    /// True if node A is an ancestor of node B, false otherwise.
    pub fn is_a_ancestor_of_b(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }

        let mut current_id = b;
        while current_id > a {
            current_id = self.nodes[current_id as usize].parent_id;
        }

        current_id == a
    }

    /// Finds the lowest common ancestor of two nodes.
    ///
    /// # Arguments
    ///
    /// * `a` - The taxonomy ID of node A.
    /// * `b` - The taxonomy ID of node B.
    ///
    /// # Returns
    ///
    /// The taxonomy ID of the lowest common ancestor of nodes A and B.
    pub fn lowest_common_ancestor(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }

        let mut current_a = a;
        let mut current_b = b;
        while current_a != current_b {
            if current_a > current_b {
                current_a = self.nodes[current_a as usize].parent_id;
            } else {
                current_b = self.nodes[current_b as usize].parent_id;
            }
        }

        current_a
    }

    /// Writes the taxonomy data to a file.
    ///
    /// # Arguments
    ///
    /// * `filename` - The path to the output file.
    pub fn write_to_disk<P: AsRef<Path>>(&self, filename: P) {
        let mut file = File::create(filename).expect("Failed to create output file");

        file.write_all(FILE_MAGIC)
            .expect("Failed to write file magic");

        let node_count = self.nodes.len() as u64;
        let name_data_len = self.name_data.len() as u64;
        let rank_data_len = self.rank_data.len() as u64;

        bincode::serialize_into(&mut file, &node_count).expect("Failed to write node count");
        bincode::serialize_into(&mut file, &name_data_len)
            .expect("Failed to write name data length");
        bincode::serialize_into(&mut file, &rank_data_len)
            .expect("Failed to write rank data length");

        bincode::serialize_into(&mut file, &self.nodes).expect("Failed to write taxonomy nodes");
        file.write_all(self.name_data.as_bytes())
            .expect("Failed to write name data");
        file.write_all(self.rank_data.as_bytes())
            .expect("Failed to write rank data");
    }

    /// Generates a mapping from external IDs to internal IDs.
    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);

        for (i, node) in self.nodes.iter().enumerate() {
            self.external_to_internal_id_map
                .insert(node.external_id, i as u64);
        }
    }

    /// Retrieves the internal ID for a given external ID.
    ///
    /// # Arguments
    ///
    /// * `external_id` - The external taxonomy ID.
    ///
    /// # Returns
    ///
    /// The corresponding internal ID, or 0 if not found.
    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        self.external_to_internal_id_map
            .get(&external_id)
            .copied()
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const NODES_DATA: &str =
        "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
3\t|\t1\t|\tclade\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|";

    const NAMES_DATA: &str = "1\t|\tall\t|\t\t|\tsynonym\t|
1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|
3\t|\tAT-rich\t|\t\t|\tsynonym\t|";

    #[test]
    fn test_ncbi_taxonomy_new() {
        let nodes_file = Cursor::new(NODES_DATA);
        let names_file = Cursor::new(NAMES_DATA);

        let taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
        println!("{:#?}", taxonomy);

        assert_eq!(taxonomy.parent_map.len(), 3);
        assert_eq!(taxonomy.name_map.len(), 2);
        assert_eq!(taxonomy.rank_map.len(), 3);
        assert_eq!(taxonomy.child_map.len(), 1);
        assert_eq!(taxonomy.marked_nodes.len(), 1);
        assert_eq!(taxonomy.known_ranks.len(), 3);

        assert_eq!(taxonomy.parent_map[&1], 0);
        assert_eq!(taxonomy.parent_map[&2], 1);
        assert_eq!(taxonomy.parent_map[&3], 1);

        assert_eq!(taxonomy.name_map[&1], "root");
        assert_eq!(taxonomy.name_map[&2], "Bacteria");
        // assert_eq!(taxonomy.name_map[&3], "AT-rich"); // Not "scientific name"

        assert_eq!(taxonomy.rank_map[&1], "no rank");
        assert_eq!(taxonomy.rank_map[&2], "superkingdom");
        assert_eq!(taxonomy.rank_map[&3], "clade");

        assert!(taxonomy.child_map[&1].contains(&2));
        assert!(taxonomy.child_map[&1].contains(&3));

        assert!(taxonomy.marked_nodes.contains(&1));

        assert!(taxonomy.known_ranks.contains("no rank"));
        assert!(taxonomy.known_ranks.contains("superkingdom"));
        assert!(taxonomy.known_ranks.contains("clade"));
    }

    #[test]
    fn test_ncbi_taxonomy_mark_node() {
        let nodes_file = Cursor::new(NODES_DATA);
        let names_file = Cursor::new(NAMES_DATA);

        let mut taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
        println!("{:?}", taxonomy);

        taxonomy.mark_node(3);

        assert!(taxonomy.marked_nodes.contains(&1));
        assert!(taxonomy.marked_nodes.contains(&3));
    }

    #[test]
    fn test_taxonomy_is_a_ancestor_of_b() {
        let nodes_file = Cursor::new(NODES_DATA);
        let names_file = Cursor::new(NAMES_DATA);

        let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
        ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

        let mut taxonomy = Taxonomy::new("test_taxonomy.k2d");
        taxonomy.generate_external_to_internal_id_map();

        let root_id = taxonomy.get_internal_id(1);
        let bacteria_id = taxonomy.get_internal_id(2);
        let at_rich_id = taxonomy.get_internal_id(3);

        assert!(taxonomy.is_a_ancestor_of_b(root_id, bacteria_id));
        assert!(taxonomy.is_a_ancestor_of_b(root_id, at_rich_id));
        assert!(!taxonomy.is_a_ancestor_of_b(bacteria_id, root_id));
        assert!(!taxonomy.is_a_ancestor_of_b(bacteria_id, at_rich_id));
        assert!(!taxonomy.is_a_ancestor_of_b(at_rich_id, root_id));
        assert!(!taxonomy.is_a_ancestor_of_b(at_rich_id, bacteria_id));
    }

    #[test]
    fn test_taxonomy_lowest_common_ancestor() {
        let nodes_file = Cursor::new(NODES_DATA);
        let names_file = Cursor::new(NAMES_DATA);

        let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
        ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

        let taxonomy = Taxonomy::new("test_taxonomy.k2d");

        assert_eq!(taxonomy.lowest_common_ancestor(1, 2), 1);
        assert_eq!(taxonomy.lowest_common_ancestor(1, 3), 1);
        assert_eq!(taxonomy.lowest_common_ancestor(2, 3), 1);
        assert_eq!(taxonomy.lowest_common_ancestor(2, 2), 2);
        assert_eq!(taxonomy.lowest_common_ancestor(3, 3), 3);
    }

    #[test]
    fn test_taxonomy_get_internal_id() {
        let nodes_file = Cursor::new(NODES_DATA);
        let names_file = Cursor::new(NAMES_DATA);

        let ncbi_taxonomy = NCBITaxonomy::from_readers(nodes_file, names_file);
        ncbi_taxonomy.convert_to_kraken_taxonomy("test_taxonomy.k2d");

        let mut taxonomy = Taxonomy::new("test_taxonomy.k2d");
        taxonomy.generate_external_to_internal_id_map();

        assert_eq!(taxonomy.get_internal_id(1), 1);
        assert_eq!(taxonomy.get_internal_id(2), 2);
        assert_eq!(taxonomy.get_internal_id(3), 3);
        assert_eq!(taxonomy.get_internal_id(4), 0);
    }
}
