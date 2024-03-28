/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use anyhow::{anyhow, Context, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
// use crate::err::Result;
// use crate::taxnode::TaxonomyNode;
// use crate::utilities::{err, errx, EX_DATAERR, EX_NOINPUT, EX_OSERR};

const FILE_MAGIC: &str = "KRAKEN_DB_2";

pub struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

impl Default for TaxonomyNode {
    fn default() -> Self {
        Self {
            parent_id: 0,
            first_child: 0,
            child_count: 0,
            name_offset: 0,
            rank_offset: 0,
            external_id: 0,
            godparent_id: 0,
        }
    }
}

pub struct NCBITaxonomy {
    parent_map: HashMap<u64, u64>,
    child_map: HashMap<u64, HashSet<u64>>,
    rank_map: HashMap<u64, String>,
    name_map: HashMap<u64, String>,
    known_ranks: HashSet<String>,
    marked_nodes: HashSet<u64>,
}

impl NCBITaxonomy {
    pub fn new(nodes_filename: &str, names_filename: &str) -> Result<Self> {
        let mut child_map: HashMap<u64, HashSet<u64>> = HashMap::new();
        for line in BufReader::new(nodes_file).lines() {
            let line = line?;
            let mut fields = line.trim_end_matches("\t|").split(delim);

            let node_id: u64 = fields
                .next()
                .context("missing node id in taxonomy file")
                .parse()
                .context("invalid node id");
            let parent_id: u64 = fields.next().unwrap_or("0").parse()?;
            let rank = fields.next().unwrap_or("").to_string();

            if node_id == 1 {
                parent_map.insert(node_id, 0);
            } else {
                parent_map.insert(node_id, parent_id);
            }
            child_map.entry(parent_id).or_default().insert(node_id);
            rank_map.insert(node_id, rank.clone());
            known_ranks.insert(rank);
        }

        for line in BufReader::new(names_file).lines() {
            let line = line?;
            let mut fields = line.trim_end_matches("\t|").split(delim);

            let node_id: u64 = fields
                .next()
                .context("attempt to create taxonomy w/ node ID == 0")
                .parse()
                .context("invalid node id");
            let name = fields.next().unwrap_or("").to_string();
            let name_type = fields.nth(1).unwrap_or("");

            if name_type == "scientific name" {
                name_map.insert(node_id, name);
            }
        }

        marked_nodes.insert(1); // mark root node

        Ok(Self {
            parent_map,
            child_map,
            rank_map,
            name_map,
            known_ranks,
            marked_nodes,
        })
    }

    // Mark the given taxonomy node and all its unmarked ancestors
    pub fn mark_node(&mut self, taxid: u64) {
        let mut taxid = taxid;
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            taxid = self.parent_map[&taxid];
        }
    }

    pub fn convert_to_kraken_taxonomy<P: AsRef<Path>>(&self, filename: P) -> Result<()> {
        let mut taxo = Taxonomy::default();

        taxo.node_count = self.marked_nodes.len() + 1; // +1 because 0 is illegal value
        taxo.nodes = vec![TaxonomyNode::default(); taxo.node_count];

        let mut name_data = String::new();
        let mut rank_data = String::new();

        // Because so many of the node rank names are shared, we only store one copy
        // of each rank
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
        while !bfs_queue.is_empty() {
            internal_node_id += 1;
            let external_node_id = bfs_queue.pop_front().unwrap();
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = TaxonomyNode::default();
            node.parent_id = external_id_map[&self.parent_map[&external_node_id]];
            node.external_id = external_node_id;
            node.rank_offset = rank_offsets[&self.rank_map[&external_node_id]];
            node.name_offset = name_data.len() as u64;
            node.first_child = internal_node_id + 1 + bfs_queue.len() as u64;

            for child_node in self
                .child_map
                .get(&external_node_id)
                .unwrap_or(&HashSet::new())
            {
                // Only add marked nodes to our internal tree
                if self.marked_nodes.contains(child_node) {
                    bfs_queue.push_back(*child_node);
                    node.child_count += 1;
                }
            }
            taxo.nodes[internal_node_id as usize] = node;

            let name = self
                .name_map
                .get(&external_node_id)
                .unwrap_or(&String::new());
            name_data.push_str(name);
            name_data.push('\0');
        }

        taxo.rank_data = rank_data.into_bytes();
        taxo.name_data = name_data.into_bytes();

        taxo.write_to_disk(filename)
    }
}

pub struct Taxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: Vec<u8>,
    pub rank_data: Vec<u8>,
    pub node_count: usize,
    pub file_backed: bool,
    pub taxonomy_data_file: Option<File>,
    pub external_to_internal_id_map: HashMap<u64, usize>,
}

impl Default for Taxonomy {
    fn default() -> Self {
        Self {
            nodes: Vec::new(),
            name_data: Vec::new(),
            rank_data: Vec::new(),
            node_count: 0,
            file_backed: false,
            taxonomy_data_file: None,
            external_to_internal_id_map: HashMap::new(),
        }
    }
}

impl Taxonomy {
    pub fn new<P: AsRef<Path>>(filename: P, memory_mapping: bool) -> Result<Self> {
        Self::init(filename.as_ref(), memory_mapping)
    }

    fn init(filename: &Path, memory_mapping: bool) -> Result<Self> {
        let mut taxo = Self::default();

        if memory_mapping {
            let file = File::open(filename)?;
            let mut reader = BufReader::new(&file);

            let mut magic = [0; 11];
            reader.read_exact(&mut magic)?;
            if magic != FILE_MAGIC.as_bytes() {
                return Err(anyhow!(
                    "attempt to load taxonomy from malformed file {}",
                    filename.display()
                )
                .context(format!(
                    "Failed to load taxonomy from file {}",
                    filename.display()
                )));
            }

            reader.read_exact(&mut taxo.node_count.to_le_bytes())?;
            let name_data_len = reader.read_u64::<LittleEndian>()?;
            let rank_data_len = reader.read_u64::<LittleEndian>()?;

            let mut nodes_data = vec![0; std::mem::size_of::<TaxonomyNode>() * taxo.node_count];
            reader.read_exact(&mut nodes_data)?;
            taxo.nodes = unsafe {
                std::slice::from_raw_parts(
                    nodes_data.as_ptr() as *const TaxonomyNode,
                    taxo.node_count,
                )
            }
            .to_vec();

            let mut name_data = vec![0; name_data_len as usize];
            reader.read_exact(&mut name_data)?;
            taxo.name_data = name_data;

            let mut rank_data = vec![0; rank_data_len as usize];
            reader.read_exact(&mut rank_data)?;
            taxo.rank_data = rank_data;

            taxo.file_backed = true;
            taxo.taxonomy_data_file = Some(file);
        } else {
            let mut file = File::open(filename)?;

            let mut magic = [0; 11];
            file.read_exact(&mut magic)?;
            if magic != FILE_MAGIC.as_bytes() {
                return Err(
                    anyhow!("malformed tax file {}", filename.display()).context(format!(
                        "Failed to load taxonomy from file {}",
                        filename.display()
                    )),
                );
            }

            file.read_exact(unsafe {
                std::slice::from_raw_parts_mut(
                    &mut taxo.node_count as *mut usize as *mut u8,
                    std::mem::size_of::<usize>(),
                )
            })?;
            let name_data_len = file.read_u64::<LittleEndian>()?;
            let rank_data_len = file.read_u64::<LittleEndian>()?;

            taxo.nodes = vec![TaxonomyNode::default(); taxo.node_count];
            unsafe {
                file.read_exact(std::slice::from_raw_parts_mut(
                    taxo.nodes.as_mut_ptr() as *mut u8,
                    std::mem::size_of::<TaxonomyNode>() * taxo.node_count,
                ))?;
            }

            taxo.name_data = vec![0; name_data_len as usize];
            file.read_exact(&mut taxo.name_data)?;

            taxo.rank_data = vec![0; rank_data_len as usize];
            file.read_exact(&mut taxo.rank_data)?;

            if file.read(&mut [0]).is_ok() {
                return Err(anyhow!(
                    "attempt to load taxonomy from malformed file {}",
                    filename.display()
                )
                .context(format!(
                    "Failed to load taxonomy from file {}",
                    filename.display()
                )));
            }
        }

        Ok(taxo)
    }

    // Logic here depends on higher nodes having smaller IDs
    // Idea: advance B tracker up tree, A is ancestor iff B tracker hits A
    pub fn is_ancestor(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        let mut b = b;
        while b > a {
            b = self.nodes[b as usize].parent_id;
        }
        b == a
    }

    // Logic here depends on higher nodes having smaller IDs
    // Idea: track two nodes, advance lower tracker up tree, trackers meet @ LCA
    pub fn lowest_common_ancestor(&self, mut a: u64, mut b: u64) -> u64 {
        if a == 0 || b == 0 {
            // LCA(x,0) = LCA(0,x) = x
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

    // Dump binary data to file
    fn write_to_disk<P: AsRef<Path>>(&self, filename: P) -> Result<()> {
        let mut file = File::create(filename)?;
        file.write_all(FILE_MAGIC.as_bytes())?;
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &self.node_count as *const usize as *const u8,
                std::mem::size_of::<usize>(),
            )
        })?;
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &(self.name_data.len() as u64) as *const u64 as *const u8,
                std::mem::size_of::<u64>(),
            )
        })?;
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                &(self.rank_data.len() as u64) as *const u64 as *const u8,
                std::mem::size_of::<u64>(),
            )
        })?;
        file.write_all(unsafe {
            std::slice::from_raw_parts(
                self.nodes.as_ptr() as *const u8,
                std::mem::size_of::<TaxonomyNode>() * self.node_count,
            )
        })?;
        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        if file.sync_all().is_err() {
            return Err(anyhow!(
                "attempt to load taxonomy from malformed file {}",
                filename.display()
            )
            .context(format!(
                "Failed to load taxonomy from file {}",
                filename.display()
            )));
        }

        Ok(())
    }

    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);
        for (i, node) in self.nodes.iter().enumerate().skip(1) {
            self.external_to_internal_id_map.insert(node.external_id, i);
        }
    }
    pub fn move_to_memory(&mut self) {
        // Implementation for moving taxonomy data to memory
        // (Not provided in the original source code)
        todo!("Implementation missing");
    }

    pub fn get_internal_id(&self, external_id: u64) -> Option<usize> {
        self.external_to_internal_id_map.get(&external_id).copied()
    }
}
