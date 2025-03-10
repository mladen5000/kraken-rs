// Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
//
// This file is part of the Kraken 2 taxonomic sequence classification system.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use crate::mmap_file::MMapFile;

#[derive(Debug, Clone, Copy, Default)]
pub struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

#[derive(Clone)]
pub struct NCBITaxonomyImpl {
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

impl NCBITaxonomyImpl {
    pub fn new<P: AsRef<Path>>(nodes_filename: P, names_filename: P) -> io::Result<Self> {
        let mut parent_map = HashMap::new();
        let mut name_map = HashMap::new();
        let mut rank_map = HashMap::new();
        let mut child_map = HashMap::new();
        let mut marked_nodes = HashSet::new();
        let mut known_ranks = HashSet::new();

        // Parse nodes file
        let nodes_file = File::open(nodes_filename)?;
        let nodes_reader = BufReader::new(nodes_file);

        for line in nodes_reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }

            let line = line.trim_end_matches("\t|");
            let delim = "\t|\t";
            let parts: Vec<&str> = line.split(delim).collect();

            if parts.len() < 3 {
                continue;
            }

            let node_id = parts[0]
                .parse::<u64>()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid node ID"))?;

            if node_id == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Attempt to create taxonomy with node ID == 0",
                ));
            }

            let mut parent_id = parts[1]
                .parse::<u64>()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid parent ID"))?;

            let rank = parts[2].to_string();

            // Root node's parent is 0
            if node_id == 1 {
                parent_id = 0;
            }

            parent_map.insert(node_id, parent_id);
            child_map
                .entry(parent_id)
                .or_insert_with(HashSet::new)
                .insert(node_id);
            rank_map.insert(node_id, rank.clone());
            known_ranks.insert(rank);
        }

        // Parse names file
        let names_file = File::open(names_filename)?;
        let names_reader = BufReader::new(names_file);

        for line in names_reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }

            let line = line.trim_end_matches("\t|");
            let delim = "\t|\t";
            let parts: Vec<&str> = line.split(delim).collect();

            if parts.len() < 4 {
                continue;
            }

            let node_id = parts[0]
                .parse::<u64>()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid node ID"))?;

            if node_id == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Node ID 0 not allowed",
                ));
            }

            let name = parts[1].to_string();
            let name_type = parts[3];

            // Only record scientific names
            if name_type == "scientific name" {
                name_map.insert(node_id, name);
            }
        }

        // Mark root node
        marked_nodes.insert(1);

        Ok(Self {
            parent_map,
            name_map,
            rank_map,
            child_map,
            marked_nodes,
            known_ranks,
        })
    }

    pub fn mark_node(&mut self, mut taxid: u64) {
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            taxid = *self.parent_map.get(&taxid).unwrap_or(&0);
        }
    }

    pub fn convert_to_kraken_taxonomy<P: AsRef<Path>>(&self, filename: P) -> io::Result<()> {
        let mut taxo = Taxonomy::default();
        let zeroes_node = TaxonomyNode::default();

        // +1 because 0 is illegal value
        taxo.node_count = self.marked_nodes.len() + 1;
        let mut nodes = vec![zeroes_node; taxo.node_count];

        let mut name_data = String::new();
        let mut rank_data = String::new();

        let mut rank_offsets = HashMap::new();
        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.push_str(rank);
            rank_data.push('\0');
        }

        let mut internal_node_id = 0;
        let mut external_id_map = HashMap::new();
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1);

        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(1);

        while !bfs_queue.is_empty() {
            internal_node_id += 1;
            let external_node_id = bfs_queue.pop_front().unwrap();
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = zeroes_node;
            node.parent_id = external_id_map[&self.parent_map[&external_node_id]];
            node.external_id = external_node_id;
            node.rank_offset = *rank_offsets
                .get(&self.rank_map[&external_node_id])
                .unwrap_or(&0);
            node.name_offset = name_data.len() as u64;
            node.first_child = internal_node_id + 1 + bfs_queue.len() as u64;

            if let Some(children) = self.child_map.get(&external_node_id) {
                for &child_node in children {
                    if self.marked_nodes.contains(&child_node) {
                        bfs_queue.push_back(child_node);
                        node.child_count += 1;
                    }
                }
            }

            nodes[internal_node_id as usize] = node;

            if let Some(name) = self.name_map.get(&external_node_id) {
                name_data.push_str(name);
                name_data.push('\0');
            }
        }

        taxo.nodes = nodes;
        taxo.rank_data = rank_data.into_bytes();
        taxo.rank_data_len = taxo.rank_data.len();
        taxo.name_data = name_data.into_bytes();
        taxo.name_data_len = taxo.name_data.len();

        taxo.write_to_disk(filename)
    }
}

#[derive(Default)]
pub struct Taxonomy {
    pub file_backed: bool,
    pub nodes: Vec<TaxonomyNode>,
    pub node_count: usize,
    pub name_data: Vec<u8>,
    pub name_data_len: usize,
    pub rank_data: Vec<u8>,
    pub rank_data_len: usize,
    pub external_to_internal_id_map: HashMap<u64, u64>,
    pub taxonomy_data_file: Option<MMapFile>,
}

impl Taxonomy {
    const FILE_MAGIC: &'static str = "K2TAXDAT";

    pub fn from_ncbi_dmp(nodes_path: &str, names_path: &str) -> io::Result<Self> {
        let ncbi_tax = NCBITaxonomyImpl::new(nodes_path, names_path)?;
        ncbi_tax.convert_to_kraken_taxonomy(nodes_path)?;
        Ok(Self::default())
    }

    pub fn from_kraken_taxonomy(path: &str) -> io::Result<Self> {
        Self::new(path, false)
    }

    pub fn new<P: AsRef<Path>>(filename: P, memory_mapping: bool) -> io::Result<Self> {
        if memory_mapping {
            Self::init_memory_mapped(filename)
        } else {
            Self::init_in_memory(filename)
        }
    }

    fn init_memory_mapped<P: AsRef<Path>>(filename: P) -> io::Result<Self> {
        let taxonomy_data_file =
            MMapFile::open_file(filename).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let ptr = taxonomy_data_file.fptr();

        // Check file magic
        let magic_len = Self::FILE_MAGIC.len();
        let magic = std::str::from_utf8(&ptr[..magic_len])
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid file format"))?;

        if magic != Self::FILE_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Malformed taxonomy file",
            ));
        }

        let mut offset = magic_len;

        // Read node count
        let node_count = unsafe {
            let ptr = ptr[offset..].as_ptr() as *const usize;
            *ptr
        };
        offset += std::mem::size_of::<usize>();

        // Read name data length
        let name_data_len = unsafe {
            let ptr = ptr[offset..].as_ptr() as *const usize;
            *ptr
        };
        offset += std::mem::size_of::<usize>();

        // Read rank data length
        let rank_data_len = unsafe {
            let ptr = ptr[offset..].as_ptr() as *const usize;
            *ptr
        };
        offset += std::mem::size_of::<usize>();

        // Get nodes pointer
        let nodes_size = std::mem::size_of::<TaxonomyNode>() * node_count;
        let nodes_data = &ptr[offset..(offset + nodes_size)];
        let nodes: Vec<TaxonomyNode> = unsafe {
            std::slice::from_raw_parts(nodes_data.as_ptr() as *const TaxonomyNode, node_count)
                .to_vec()
        };
        offset += nodes_size;

        // Get name data
        let name_data = ptr[offset..(offset + name_data_len)].to_vec();
        offset += name_data_len;

        // Get rank data
        let rank_data = ptr[offset..(offset + rank_data_len)].to_vec();

        Ok(Self {
            file_backed: true,
            nodes,
            node_count,
            name_data,
            name_data_len,
            rank_data,
            rank_data_len,
            external_to_internal_id_map: HashMap::new(),
            taxonomy_data_file: Some(taxonomy_data_file),
        })
    }

    fn init_in_memory<P: AsRef<Path>>(filename: P) -> io::Result<Self> {
        let mut file = File::open(filename)?;

        // Check file magic
        let magic_len = Self::FILE_MAGIC.len();
        let mut magic = vec![0; magic_len];
        file.read_exact(&mut magic)?;

        let magic_str = std::str::from_utf8(&magic)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid file format"))?;

        if magic_str != Self::FILE_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Malformed taxonomy file",
            ));
        }

        // Read node count
        let mut node_count_bytes = [0; std::mem::size_of::<usize>()];
        file.read_exact(&mut node_count_bytes)?;
        let node_count = usize::from_ne_bytes(node_count_bytes);

        // Read name data length
        let mut name_data_len_bytes = [0; std::mem::size_of::<usize>()];
        file.read_exact(&mut name_data_len_bytes)?;
        let name_data_len = usize::from_ne_bytes(name_data_len_bytes);

        // Read rank data length
        let mut rank_data_len_bytes = [0; std::mem::size_of::<usize>()];
        file.read_exact(&mut rank_data_len_bytes)?;
        let rank_data_len = usize::from_ne_bytes(rank_data_len_bytes);

        // Read nodes
        let mut nodes = vec![TaxonomyNode::default(); node_count];
        for i in 0..node_count {
            // Read each field of the node
            let mut parent_id_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut parent_id_bytes)?;
            nodes[i].parent_id = u64::from_ne_bytes(parent_id_bytes);

            let mut first_child_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut first_child_bytes)?;
            nodes[i].first_child = u64::from_ne_bytes(first_child_bytes);

            let mut child_count_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut child_count_bytes)?;
            nodes[i].child_count = u64::from_ne_bytes(child_count_bytes);

            let mut name_offset_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut name_offset_bytes)?;
            nodes[i].name_offset = u64::from_ne_bytes(name_offset_bytes);

            let mut rank_offset_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut rank_offset_bytes)?;
            nodes[i].rank_offset = u64::from_ne_bytes(rank_offset_bytes);

            let mut external_id_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut external_id_bytes)?;
            nodes[i].external_id = u64::from_ne_bytes(external_id_bytes);

            let mut godparent_id_bytes = [0; std::mem::size_of::<u64>()];
            file.read_exact(&mut godparent_id_bytes)?;
            nodes[i].godparent_id = u64::from_ne_bytes(godparent_id_bytes);
        }

        // Read name data
        let mut name_data = vec![0; name_data_len];
        file.read_exact(&mut name_data)?;

        // Read rank data
        let mut rank_data = vec![0; rank_data_len];
        file.read_exact(&mut rank_data)?;

        if file.read(&mut [0])? != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Extra data in taxonomy file",
            ));
        }

        Ok(Self {
            file_backed: false,
            nodes,
            node_count,
            name_data,
            name_data_len,
            rank_data,
            rank_data_len,
            external_to_internal_id_map: HashMap::new(),
            taxonomy_data_file: None,
        })
    }

    pub fn write_to_disk<P: AsRef<Path>>(&self, filename: P) -> io::Result<()> {
        let mut file = File::create(filename)?;

        // Write magic identifier
        file.write_all(Self::FILE_MAGIC.as_bytes())?;

        // Write node count
        file.write_all(&self.node_count.to_ne_bytes())?;

        // Write name data length
        file.write_all(&self.name_data_len.to_ne_bytes())?;

        // Write rank data length
        file.write_all(&self.rank_data_len.to_ne_bytes())?;

        // Write nodes
        for node in &self.nodes {
            file.write_all(&node.parent_id.to_ne_bytes())?;
            file.write_all(&node.first_child.to_ne_bytes())?;
            file.write_all(&node.child_count.to_ne_bytes())?;
            file.write_all(&node.name_offset.to_ne_bytes())?;
            file.write_all(&node.rank_offset.to_ne_bytes())?;
            file.write_all(&node.external_id.to_ne_bytes())?;
            file.write_all(&node.godparent_id.to_ne_bytes())?;
        }

        // Write name and rank data
        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        Ok(())
    }

    pub fn nodes(&self) -> &[TaxonomyNode] {
        &self.nodes
    }

    pub fn name_data(&self) -> &[u8] {
        &self.name_data
    }

    pub fn rank_data(&self) -> &[u8] {
        &self.rank_data
    }

    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        *self
            .external_to_internal_id_map
            .get(&external_id)
            .unwrap_or(&0)
    }

    pub fn lowest_common_ancestor(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            // LCA(x,0) = LCA(0,x) = x
            return if a == 0 { b } else { a };
        }
        let mut a = a;
        let mut b = b;
        while a != b {
            if a > b {
                a = self.nodes[a as usize].parent_id;
            } else {
                b = self.nodes[b as usize].parent_id;
            }
        }
        a
    }

    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);
        for i in 1..self.node_count {
            self.external_to_internal_id_map
                .insert(self.nodes[i].external_id, i as u64);
        }
    }

    pub fn mark_node(&self, _taxid: u64) {
        // Only needed in NCBITaxonomy
    }

    pub fn is_a_ancestor_of_b(&self, potential_ancestor: u64, potential_descendant: u64) -> bool {
        if potential_ancestor == 0 || potential_descendant == 0 {
            return false;
        }
        // Logic here depends on higher nodes having smaller IDs
        // Advance descendant tracker up tree, potential_ancestor is ancestor iff tracker hits it
        let mut b = potential_descendant;
        while b > potential_ancestor {
            b = self.nodes[b as usize].parent_id;
        }
        b == potential_ancestor
    }

    pub fn name_at(&self, offset: u64) -> Option<&str> {
        // Find the end of the name string (marked by null terminator)
        let start = offset as usize;
        if start >= self.name_data.len() {
            return None;
        }
        let mut end = start;
        while end < self.name_data.len() && self.name_data[end] != 0 {
            end += 1;
        }
        std::str::from_utf8(&self.name_data[start..end]).ok()
    }

    pub fn rank_data_at(&self, offset: u64) -> &str {
        // Find the end of the rank string (marked by null terminator)
        let start = offset as usize;
        if start >= self.rank_data.len() {
            return "";
        }
        let mut end = start;
        while end < self.rank_data.len() && self.rank_data[end] != 0 {
            end += 1;
        }
        std::str::from_utf8(&self.rank_data[start..end]).unwrap_or("")
    }

    pub fn name_data_at(&self, offset: u64) -> &str {
        // Find the end of the name string (marked by null terminator)
        let start = offset as usize;
        if start >= self.name_data.len() {
            return "";
        }
        let mut end = start;
        while end < self.name_data.len() && self.name_data[end] != 0 {
            end += 1;
        }
        std::str::from_utf8(&self.name_data[start..end]).unwrap_or("")
    }
}

impl Drop for Taxonomy {
    fn drop(&mut self) {
        self.taxonomy_data_file.take(); // Drop the mmap file if it exists
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{BufWriter, Write};
    use tempfile::tempdir;

    fn create_test_taxonomy() -> (Taxonomy, tempfile::TempDir) {
        let temp_dir = tempdir().unwrap();
        let taxonomy_path = temp_dir.path().join("taxonomy.bin");

        // Create a simple taxonomy with known structure:
        //       1 (root)
        //      / \
        //     2   3
        //    /     \
        //   4       5
        let mut taxonomy = Taxonomy::default();
        taxonomy.node_count = 6; // 0-5, where 0 is unused

        let mut nodes = vec![TaxonomyNode::default(); 6];

        // Root node (id 1)
        nodes[1] = TaxonomyNode {
            parent_id: 0,
            first_child: 2,
            child_count: 2,
            name_offset: 0,
            rank_offset: 0,
            external_id: 1,
            godparent_id: 0,
        };

        // Node 2
        nodes[2] = TaxonomyNode {
            parent_id: 1,
            first_child: 4,
            child_count: 1,
            name_offset: 5,
            rank_offset: 0,
            external_id: 2,
            godparent_id: 0,
        };

        // Node 3
        nodes[3] = TaxonomyNode {
            parent_id: 1,
            first_child: 5,
            child_count: 1,
            name_offset: 10,
            rank_offset: 0,
            external_id: 3,
            godparent_id: 0,
        };

        // Node 4
        nodes[4] = TaxonomyNode {
            parent_id: 2,
            first_child: 0,
            child_count: 0,
            name_offset: 15,
            rank_offset: 0,
            external_id: 4,
            godparent_id: 0,
        };

        // Node 5
        nodes[5] = TaxonomyNode {
            parent_id: 3,
            first_child: 0,
            child_count: 0,
            name_offset: 20,
            rank_offset: 0,
            external_id: 5,
            godparent_id: 0,
        };

        taxonomy.nodes = nodes;

        // Create name data
        taxonomy.name_data = b"root\0node2\0node3\0node4\0node5\0".to_vec();
        taxonomy.name_data_len = taxonomy.name_data.len();

        // Create rank data
        taxonomy.rank_data = b"root\0phylum\0class\0genus\0species\0".to_vec();
        taxonomy.rank_data_len = taxonomy.rank_data.len();

        // Create external to internal ID map
        taxonomy.external_to_internal_id_map.insert(0, 0);
        taxonomy.external_to_internal_id_map.insert(1, 1);
        taxonomy.external_to_internal_id_map.insert(2, 2);
        taxonomy.external_to_internal_id_map.insert(3, 3);
        taxonomy.external_to_internal_id_map.insert(4, 4);
        taxonomy.external_to_internal_id_map.insert(5, 5);

        // Write the taxonomy to file for tests that need it
        taxonomy.write_to_disk(&taxonomy_path).unwrap();

        (taxonomy, temp_dir)
    }

    #[test]
    fn test_lowest_common_ancestor() {
        let (taxonomy, _temp_dir) = create_test_taxonomy();

        // Test with node being its own ancestor
        assert_eq!(taxonomy.lowest_common_ancestor(2, 2), 2);

        // Test with one node being 0
        assert_eq!(taxonomy.lowest_common_ancestor(0, 3), 3);
        assert_eq!(taxonomy.lowest_common_ancestor(4, 0), 4);

        // Test with siblings
        assert_eq!(taxonomy.lowest_common_ancestor(2, 3), 1);

        // Test with one node being ancestor of the other
        assert_eq!(taxonomy.lowest_common_ancestor(1, 4), 1);
        assert_eq!(taxonomy.lowest_common_ancestor(4, 2), 2);

        // Test with nodes in different subtrees
        assert_eq!(taxonomy.lowest_common_ancestor(4, 5), 1);
    }

    #[test]
    fn test_is_a_ancestor_of_b() {
        let (taxonomy, _temp_dir) = create_test_taxonomy();

        // Test direct parent-child relationship
        assert!(taxonomy.is_a_ancestor_of_b(2, 4));
        assert!(taxonomy.is_a_ancestor_of_b(1, 4));

        // Test non-ancestor relationship
        assert!(!taxonomy.is_a_ancestor_of_b(2, 5));
        assert!(!taxonomy.is_a_ancestor_of_b(4, 2)); // Child is not ancestor of parent

        // Test self-relationship
        assert!(taxonomy.is_a_ancestor_of_b(3, 3));

        // Test with zero values
        assert!(!taxonomy.is_a_ancestor_of_b(0, 1));
        assert!(!taxonomy.is_a_ancestor_of_b(1, 0));
    }

    #[test]
    fn test_name_at() {
        let (taxonomy, _temp_dir) = create_test_taxonomy();

        // Test valid offsets
        assert_eq!(taxonomy.name_at(0), Some("root"));
        assert_eq!(taxonomy.name_at(5), Some("node2"));
        assert_eq!(taxonomy.name_at(15), Some("node4"));

        // Test invalid offset
        assert_eq!(taxonomy.name_at(100), None);
    }

    #[test]
    fn test_read_write_taxonomy() {
        let (original_taxonomy, temp_dir) = create_test_taxonomy();
        let taxonomy_path = temp_dir.path().join("test_rw.bin");

        // Write the taxonomy to a file
        original_taxonomy.write_to_disk(&taxonomy_path).unwrap();

        // Read it back
        let loaded_taxonomy = Taxonomy::new(&taxonomy_path, false).unwrap();

        // Check that the read taxonomy matches the original
        assert_eq!(loaded_taxonomy.node_count(), original_taxonomy.node_count());
        assert_eq!(
            loaded_taxonomy.name_data().len(),
            original_taxonomy.name_data().len()
        );
        assert_eq!(
            loaded_taxonomy.rank_data().len(),
            original_taxonomy.rank_data().len()
        );

        // Compare the first few nodes
        for i in 1..6 {
            assert_eq!(
                loaded_taxonomy.nodes()[i].parent_id,
                original_taxonomy.nodes()[i].parent_id
            );
            assert_eq!(
                loaded_taxonomy.nodes()[i].external_id,
                original_taxonomy.nodes()[i].external_id
            );
        }
    }

    #[test]
    fn test_ncbi_taxonomy_impl() {
        let temp_dir = tempdir().unwrap();
        let nodes_path = temp_dir.path().join("nodes.dmp");
        let names_path = temp_dir.path().join("names.dmp");

        // Create simplified NCBI taxonomy files
        {
            let mut nodes_file = BufWriter::new(File::create(&nodes_path).unwrap());
            // Format: taxid | parent taxid | rank | ...
            writeln!(nodes_file, "1\t|\t1\t|\tno rank\t|\t").unwrap(); // root
            writeln!(nodes_file, "2\t|\t1\t|\tphylum\t|\t").unwrap();
            writeln!(nodes_file, "3\t|\t1\t|\tphylum\t|\t").unwrap();
            writeln!(nodes_file, "4\t|\t2\t|\tclass\t|\t").unwrap();
            writeln!(nodes_file, "5\t|\t3\t|\tclass\t|\t").unwrap();
        }

        {
            let mut names_file = BufWriter::new(File::create(&names_path).unwrap());
            // Format: taxid | name | unique name | name class |
            writeln!(names_file, "1\t|\troot\t|\t\t|\tscientific name\t|\t").unwrap();
            writeln!(names_file, "2\t|\tBacteria\t|\t\t|\tscientific name\t|\t").unwrap();
            writeln!(names_file, "3\t|\tArchaea\t|\t\t|\tscientific name\t|\t").unwrap();
            writeln!(
                names_file,
                "4\t|\tProteobacteria\t|\t\t|\tscientific name\t|\t"
            )
            .unwrap();
            writeln!(
                names_file,
                "5\t|\tEuryarchaeota\t|\t\t|\tscientific name\t|\t"
            )
            .unwrap();
        }

        // Parse the taxonomy
        let ncbi_taxonomy = NCBITaxonomyImpl::new(&nodes_path, &names_path).unwrap();

        // Check that parsing worked
        assert_eq!(ncbi_taxonomy.parent_map.len(), 5);
        assert_eq!(ncbi_taxonomy.name_map.len(), 5);
        assert_eq!(ncbi_taxonomy.rank_map.len(), 5);

        // Test mark_node
        let mut ncbi_taxonomy_for_marking = ncbi_taxonomy.clone();
        ncbi_taxonomy_for_marking.mark_node(4);

        // This should have marked 4, 2, and 1 (the path to root)
        assert!(ncbi_taxonomy_for_marking.marked_nodes.contains(&4));
        assert!(ncbi_taxonomy_for_marking.marked_nodes.contains(&2));
        assert!(ncbi_taxonomy_for_marking.marked_nodes.contains(&1));
        assert!(!ncbi_taxonomy_for_marking.marked_nodes.contains(&3));
        assert!(!ncbi_taxonomy_for_marking.marked_nodes.contains(&5));

        // Test conversion to Kraken taxonomy
        let kraken_path = temp_dir.path().join("kraken.tax");
        ncbi_taxonomy
            .convert_to_kraken_taxonomy(&kraken_path)
            .unwrap();

        // Load the Kraken taxonomy and check it
        let kraken_taxonomy = Taxonomy::new(&kraken_path, false).unwrap();

        // We've only marked the root node in the original NCBITaxonomyImpl,
        // so only that should be in the Kraken taxonomy
        assert_eq!(kraken_taxonomy.node_count(), 2); // Node 0 + root
    }
}

pub trait NCBITaxonomyOps {
    fn from_ncbi_dmp(nodes_file: &str, names_file: &str) -> Result<Self, anyhow::Error>
    where
        Self: Sized;
    fn mark_node(&self, taxid: u64);
    fn convert_to_kraken_taxonomy(&self, output_file: &str) -> Result<(), anyhow::Error>;
}

impl NCBITaxonomyOps for Taxonomy {
    fn from_ncbi_dmp(nodes_file: &str, names_file: &str) -> Result<Self, anyhow::Error> {
        // Create new taxonomy
        let mut taxonomy = Taxonomy {
            file_backed: false,
            nodes: Vec::new(),
            node_count: 0,
            name_data: Vec::new(),
            name_data_len: 0,
            rank_data: Vec::new(),
            rank_data_len: 0,
            external_to_internal_id_map: HashMap::new(),
            taxonomy_data_file: None,
        };

        // Read nodes.dmp
        let nodes_reader = BufReader::new(File::open(nodes_file)?);
        let mut nodes_map = HashMap::new();

        for line in nodes_reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
            if fields.len() >= 3 {
                let taxid = fields[0].parse::<u64>()?;
                let parent_id = fields[1].parse::<u64>()?;
                let rank = fields[2].to_string();
                nodes_map.insert(taxid, (parent_id, rank));
            }
        }

        // Read names.dmp
        let names_reader = BufReader::new(File::open(names_file)?);
        let mut names_map = HashMap::new();

        for line in names_reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
            if fields.len() >= 4 && fields[3] == "scientific name" {
                let taxid = fields[0].parse::<u64>()?;
                let name = fields[1].to_string();
                names_map.insert(taxid, name);
            }
        }

        // Build nodes vector in order of taxids
        let mut taxids: Vec<u64> = nodes_map.keys().cloned().collect();
        taxids.sort_unstable();

        // Create mapping from external to internal IDs
        for (idx, &taxid) in taxids.iter().enumerate() {
            taxonomy
                .external_to_internal_id_map
                .insert(taxid, idx as u64);
        }

        // Build nodes vector
        for taxid in taxids {
            let (parent_id, rank) = nodes_map[&taxid].clone();
            let parent_idx = if parent_id == taxid {
                0 // Root node points to itself
            } else {
                *taxonomy
                    .external_to_internal_id_map
                    .get(&parent_id)
                    .unwrap_or(&0) as usize
            };

            let name = names_map.get(&taxid).cloned().unwrap_or_default();
            let rank_offset = taxonomy.rank_data.len() as u64;
            taxonomy.rank_data.extend_from_slice(rank.as_bytes());
            taxonomy.rank_data.push(0);
            taxonomy.rank_data_len = taxonomy.rank_data.len();

            let name_offset = taxonomy.name_data.len() as u64;
            taxonomy.name_data.extend_from_slice(name.as_bytes());
            taxonomy.name_data.push(0);
            taxonomy.name_data_len = taxonomy.name_data.len();

            taxonomy.nodes.push(TaxonomyNode {
                parent_id: parent_idx as u64,
                first_child: 0,
                child_count: 0,
                name_offset,
                rank_offset,
                external_id: taxid,
                godparent_id: 0,
            });
        }

        // Update first_child and child_count fields
        for node_idx in 0..taxonomy.nodes.len() {
            let parent_id = taxonomy.nodes[node_idx].parent_id;
            if parent_id as usize != node_idx {
                // Skip root node
                let parent = &mut taxonomy.nodes[parent_id as usize];
                if parent.first_child == 0 {
                    parent.first_child = node_idx as u64;
                }
                parent.child_count += 1;
            }
        }

        taxonomy.node_count = taxonomy.nodes.len();
        Ok(taxonomy)
    }

    fn mark_node(&self, _taxid: u64) {
        // This is a no-op in the Rust version since we don't need to mark nodes
        // during taxonomy building - we build deterministically
    }

    fn convert_to_kraken_taxonomy(&self, output_file: &str) -> Result<(), anyhow::Error> {
        let mut file = File::create(output_file)?;

        // Write magic
        file.write_all(Taxonomy::FILE_MAGIC.as_bytes())?;

        // Write header
        file.write_all(&(self.node_count as u64).to_le_bytes())?;
        file.write_all(&(self.name_data_len as u64).to_le_bytes())?;
        file.write_all(&(self.rank_data_len as u64).to_le_bytes())?;

        // Write nodes
        for node in &self.nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            file.write_all(&node.first_child.to_le_bytes())?;
            file.write_all(&node.child_count.to_le_bytes())?;
            file.write_all(&node.name_offset.to_le_bytes())?;
            file.write_all(&node.rank_offset.to_le_bytes())?;
            file.write_all(&node.external_id.to_le_bytes())?;
            file.write_all(&node.godparent_id.to_le_bytes())?;
        }

        // Write name data
        file.write_all(&self.name_data)?;

        // Write rank data
        file.write_all(&self.rank_data)?;

        Ok(())
    }
}
