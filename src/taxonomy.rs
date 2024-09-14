/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::process;
use std::str;
use std::string::String;

use serde::Serialize;

// use crate::mmap_file::MMapFile;
use mmap::{MapOption, MemoryMap};

pub type TaxId = u64;

#[derive(Debug, Default, Clone, Serialize)]
pub struct TaxonomyNode {
    pub parent_id: TaxId,    // Must be lower-numbered node
    pub first_child: TaxId,  // Must be higher-numbered node
    pub child_count: TaxId,  // Children of a node are in contiguous block
    pub name_offset: usize,  // Location of name in name data super-string
    pub rank_offset: usize,  // Location of rank in rank data super-string
    pub external_id: TaxId,  // Taxonomy ID for reporting purposes (usually NCBI)
    pub godparent_id: TaxId, // Reserved for future use to enable faster traversal
}

pub struct NCBITaxonomy {
    parent_map: HashMap<TaxId, TaxId>,
    name_map: HashMap<TaxId, String>,
    rank_map: HashMap<TaxId, String>,
    child_map: HashMap<TaxId, HashSet<TaxId>>,
    marked_nodes: HashSet<TaxId>,
    known_ranks: HashSet<String>,
}

impl NCBITaxonomy {
    pub fn new(nodes_filename: &str, names_filename: &str) -> io::Result<Self> {
        let nodes_file = File::open(nodes_filename)?;
        let names_file = File::open(names_filename)?;

        let mut parent_map = HashMap::new();
        let mut name_map = HashMap::new();
        let mut rank_map = HashMap::new();
        let mut child_map = HashMap::new();
        let mut known_ranks = HashSet::new();
        let mut marked_nodes = HashSet::new();

        // Parse nodes.dmp
        let reader = BufReader::new(nodes_file);
        for line in reader.lines() {
            let line = line?;
            let tokens: Vec<&str> = line.trim_end_matches("\t|").split("\t|\t").collect();

            if tokens.len() < 3 {
                continue;
            }

            let node_id = tokens[0].trim().parse::<TaxId>().unwrap();
            let parent_id = tokens[1].trim().parse::<TaxId>().unwrap();
            let rank = tokens[2].trim().to_string();

            if node_id == 0 {
                eprintln!("Attempt to create taxonomy with node ID == 0");
                process::exit(1);
            }

            let parent_id = if node_id == 1 { 0 } else { parent_id };

            parent_map.insert(node_id, parent_id);
            child_map
                .entry(parent_id)
                .or_insert_with(HashSet::new)
                .insert(node_id);
            rank_map.insert(node_id, rank.clone());
            known_ranks.insert(rank);
        }

        // Parse names.dmp
        let reader = BufReader::new(names_file);
        for line in reader.lines() {
            let line = line?;
            let tokens: Vec<&str> = line.trim_end_matches("\t|").split("\t|\t").collect();

            if tokens.len() < 4 {
                continue;
            }

            let node_id = tokens[0].trim().parse::<TaxId>().unwrap();
            let name = tokens[1].trim().to_string();
            let name_type = tokens[3].trim();

            if node_id == 0 {
                eprintln!("Attempt to create taxonomy with node ID == 0");
                process::exit(1);
            }

            if name_type == "scientific name" {
                name_map.insert(node_id, name);
            }
        }

        marked_nodes.insert(1); // Mark root node

        Ok(NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            marked_nodes,
            known_ranks,
        })
    }

    // Mark the given taxonomy node and all its unmarked ancestors
    pub fn mark_node(&mut self, mut taxid: TaxId) {
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            if let Some(&parent_id) = self.parent_map.get(&taxid) {
                taxid = parent_id;
            } else {
                break;
            }
        }
    }

    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> io::Result<()> {
        let mut taxo = Taxonomy::default();
        let zero_node = TaxonomyNode::default();

        taxo.node_count = self.marked_nodes.len() + 1; // +1 because 0 is illegal value
        taxo.nodes = vec![zero_node.clone(); taxo.node_count];

        let mut name_data = String::new();
        let mut rank_data = String::new();

        // Map to store rank offsets
        let mut rank_offsets = HashMap::new();
        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len());
            rank_data.push_str(rank);
            rank_data.push('\0');
        }

        let mut internal_node_id = 0;
        let mut external_id_map = HashMap::new(); // keys: ext. ID, values: int. ID
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1); // 1 is root in both NCBI and Kraken taxonomies

        // Breadth-first search through NCBI taxonomy, assigning internal IDs
        // in sequential order as nodes are encountered via BFS.
        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(1);
        while let Some(external_node_id) = bfs_queue.pop_front() {
            internal_node_id += 1;
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = zero_node.clone();
            node.parent_id = *external_id_map
                .get(&self.parent_map[&external_node_id])
                .unwrap();
            node.external_id = external_node_id;
            node.rank_offset = *rank_offsets.get(&self.rank_map[&external_node_id]).unwrap();
            node.name_offset = name_data.len();
            node.first_child = internal_node_id + bfs_queue.len() as TaxId + 1;

            if let Some(children) = self.child_map.get(&external_node_id) {
                for &child_node in children {
                    if self.marked_nodes.contains(&child_node) {
                        bfs_queue.push_back(child_node);
                        node.child_count += 1;
                    }
                }
            }

            taxo.nodes[internal_node_id as usize] = node;

            let name = self
                .name_map
                .get(&external_node_id)
                .cloned()
                .unwrap_or_else(|| String::new());
            name_data.push_str(&name);
            name_data.push('\0');
        }

        taxo.rank_data = rank_data.into_bytes();
        taxo.rank_data_len = taxo.rank_data.len();

        taxo.name_data = name_data.into_bytes();
        taxo.name_data_len = taxo.name_data.len();

        taxo.write_to_disk(filename)
    }
}

#[derive(Default)]
pub struct Taxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub node_count: usize,
    pub name_data: Vec<u8>,
    pub name_data_len: usize,
    pub rank_data: Vec<u8>,
    pub rank_data_len: usize,
    pub external_to_internal_id_map: HashMap<TaxId, TaxId>,
    file_backed: bool,
    file_magic: [u8; 8],
    taxonomy_data_file: Option<MemoryMap>,
}

impl Taxonomy {
    pub fn new(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        let mut taxo = Taxonomy {
            file_magic: *b"K2TAXDAT",
            ..Default::default()
        };
        taxo.init(filename, memory_mapping)?;
        Ok(taxo)
    }

    fn init(&mut self, filename: &str, memory_mapping: bool) -> io::Result<()> {
        if memory_mapping {
            // Open the file in read-only mode and determine its length
            let file = File::open(filename)?;
            let file_len = file.metadata()?.len() as usize;

            // Create a memory map with the length of the file
            let mmap = MemoryMap::new(file_len, &[MapOption::MapReadable])
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?; // Convert MapError to io::Error
            self.file_backed = true;
            self.taxonomy_data_file = Some(mmap); // Fixing the field name

            // Get the data from the memory-mapped file
            let data = self.taxonomy_data_file.as_ref().unwrap().data();
            let slice = unsafe { std::slice::from_raw_parts(data, file_len) };
            let mut offset = 0;

            // Validate file magic
            if &slice[offset..offset + self.file_magic.len()] != self.file_magic {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Malformed file: {}", filename),
                ));
            }
            offset += self.file_magic.len();

            // Helper to read u64 in little endian
            let read_u64 = |data: &[u8], start: usize| -> u64 {
                u64::from_le_bytes(data[start..start + 8].try_into().unwrap())
            };

            self.node_count = read_u64(slice, offset) as usize;
            offset += 8;
            self.name_data_len = read_u64(slice, offset) as usize;
            offset += 8;
            self.rank_data_len = read_u64(slice, offset) as usize;
            offset += 8;

            // Read nodes
            let node_size = std::mem::size_of::<TaxonomyNode>();
            let nodes_end = offset + self.node_count * node_size;
            if nodes_end > slice.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of file while reading nodes",
                ));
            }
            let node_slice = &slice[offset..nodes_end];
            self.nodes = unsafe {
                std::slice::from_raw_parts(
                    node_slice.as_ptr() as *const TaxonomyNode,
                    self.node_count,
                )
                .to_vec()
            };
            offset = nodes_end;

            // Read name_data
            let name_data_end = offset + self.name_data_len;
            if name_data_end > slice.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of file while reading name_data",
                ));
            }
            self.name_data = slice[offset..name_data_end].to_vec();
            offset = name_data_end;

            // Read rank_data
            let rank_data_end = offset + self.rank_data_len;
            if rank_data_end > slice.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of file while reading rank_data",
                ));
            }
            self.rank_data = slice[offset..rank_data_end].to_vec();
        } else {
            let mut file = File::open(filename)?;
            let mut buffer = Vec::new();
            file.read_to_end(&mut buffer)?;

            let mut offset = 0;

            if &buffer[offset..offset + self.file_magic.len()] != self.file_magic {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Malformed file: {}", filename),
                ));
            }
            offset += self.file_magic.len();

            self.node_count =
                u64::from_le_bytes(buffer[offset..offset + 8].try_into().unwrap()) as usize;
            offset += 8;
            self.name_data_len =
                u64::from_le_bytes(buffer[offset..offset + 8].try_into().unwrap()) as usize;
            offset += 8;
            self.rank_data_len =
                u64::from_le_bytes(buffer[offset..offset + 8].try_into().unwrap()) as usize;
            offset += 8;

            let node_slice =
                &buffer[offset..offset + self.node_count * std::mem::size_of::<TaxonomyNode>()];
            self.nodes = unsafe {
                std::slice::from_raw_parts(
                    node_slice.as_ptr() as *const TaxonomyNode,
                    self.node_count,
                )
                .to_vec()
            };
            offset += self.node_count * std::mem::size_of::<TaxonomyNode>();

            self.name_data = buffer[offset..offset + self.name_data_len].to_vec();
            offset += self.name_data_len;
            self.rank_data = buffer[offset..offset + self.rank_data_len].to_vec();
        }
        Ok(())
    }

    pub fn is_a_ancestor_of_b(&self, a: TaxId, mut b: TaxId) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        while b > a {
            b = self.nodes[b as usize].parent_id;
        }
        b == a
    }

    pub fn lowest_common_ancestor(&self, mut a: TaxId, mut b: TaxId) -> TaxId {
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }
        while a != b {
            if a > b {
                a = self.nodes[a as usize].parent_id;
            } else {
                b = self.nodes[b as usize].parent_id;
            }
        }
        a
    }

    pub fn write_to_disk(&self, filename: &str) -> io::Result<()> {
        let mut file = BufWriter::new(File::create(filename)?);

        file.write_all(&self.file_magic)?;
        file.write_all(&(self.node_count as u64).to_le_bytes())?;
        file.write_all(&(self.name_data_len as u64).to_le_bytes())?;
        file.write_all(&(self.rank_data_len as u64).to_le_bytes())?;

        for node in &self.nodes {
            let serialized_node =
                bincode::serialize(node).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            file.write_all(&serialized_node)?;
        }

        // let serialized_nodes = bincode::serialize(&self.nodes)?;
        // file.write_all(&serialized_nodes)?;

        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        Ok(())
    }

    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);
        for (i, node) in self.nodes.iter().enumerate().skip(1) {
            self.external_to_internal_id_map
                .insert(node.external_id, i as TaxId);
        }
    }

    pub fn get_internal_id(&self, external_id: TaxId) -> TaxId {
        *self
            .external_to_internal_id_map
            .get(&external_id)
            .unwrap_or(&0)
    }
}
