use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::str;

const FILE_MAGIC: &str = "TAXONOMY\n";

#[repr(C)]
#[derive(Default, Debug, Clone)]
pub struct TaxonomyNode {
    pub parent_id: u64,
    pub external_id: u64,
    pub rank_offset: u64,
    pub name_offset: u64,
    pub child_count: u32,
    pub first_child: u64,
    // Padding in C++ or no, we just keep the same fields
    // The code includes no other fields after these in snippet
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
    pub fn new(nodes_filename: &str, names_filename: &str) -> Result<NCBITaxonomy> {
        let mut parent_map = HashMap::new();
        let mut child_map: HashMap<u64, HashSet<u64>> = HashMap::new();
        let mut rank_map = HashMap::new();
        let mut known_ranks = HashSet::new();
        let mut name_map = HashMap::new();
        let mut marked_nodes = HashSet::new();

        let nodes_file = File::open(nodes_filename)
            .with_context(|| format!("error opening {}", nodes_filename))?;
        let names_file = File::open(names_filename)
            .with_context(|| format!("error opening {}", names_filename))?;

        // Parse nodes
        {
            let delim = "\t|\t";
            let reader = BufReader::new(nodes_file);
            for line_res in reader.lines() {
                let mut line = line_res?;
                // discard trailing "\t|"
                if line.ends_with("\t|") {
                    line.truncate(line.len() - 2);
                }

                let mut pos1 = 0usize;
                let mut field_ct = 0;
                let mut finished = false;
                let mut node_id: u64 = 0;
                let mut parent_id: u64 = 0;
                let mut rank = String::new();

                while field_ct < 10 && !finished {
                    field_ct += 1;
                    let pos2 = match line[pos1..].find(delim) {
                        Some(p) => p + pos1,
                        None => {
                            finished = true;
                            line.len()
                        }
                    };
                    let token = &line[pos1..pos2];
                    pos1 = if pos2 < line.len() {
                        pos2 + delim.len()
                    } else {
                        pos2
                    };

                    match field_ct {
                        1 => {
                            node_id = token
                                .parse()
                                .map_err(|_| anyhow::anyhow!("invalid node_id: {}", token))?;
                            if node_id == 0 {
                                bail!("attempt to create taxonomy w/ node ID == 0");
                            }
                        }
                        2 => {
                            parent_id = token
                                .parse()
                                .map_err(|_| anyhow::anyhow!("invalid parent_id: {}", token))?;
                        }
                        3 => {
                            rank = token.to_string();
                            finished = true;
                        }
                        _ => {}
                    }
                }

                if node_id == 1 {
                    parent_id = 0;
                }
                parent_map.insert(node_id, parent_id);
                child_map.entry(parent_id).or_default().insert(node_id);
                rank_map.insert(node_id, rank.clone());
                known_ranks.insert(rank);
            }
        }

        // Parse names
        {
            let delim = "\t|\t";
            let reader = BufReader::new(names_file);
            for line_res in reader.lines() {
                let mut line = line_res?;
                // discard trailing "\t|"
                if line.ends_with("\t|") {
                    line.truncate(line.len() - 2);
                }

                let mut pos1 = 0usize;
                let mut field_ct = 0;
                let mut finished = false;
                let mut node_id: u64 = 0;
                let mut name = String::new();

                while field_ct < 10 && !finished {
                    field_ct += 1;
                    let pos2 = match line[pos1..].find(delim) {
                        Some(p) => p + pos1,
                        None => {
                            finished = true;
                            line.len()
                        }
                    };
                    let token = &line[pos1..pos2];
                    pos1 = if pos2 < line.len() {
                        pos2 + delim.len()
                    } else {
                        pos2
                    };

                    match field_ct {
                        1 => {
                            node_id = token
                                .parse()
                                .map_err(|_| anyhow::anyhow!("invalid node_id: {}", token))?;
                            if node_id == 0 {
                                bail!("attempt to create taxonomy w/ node ID == 0");
                            }
                        }
                        2 => {
                            name = token.to_string();
                        }
                        4 => {
                            if token == "scientific name" {
                                name_map.insert(node_id, name.clone());
                            }
                            finished = true;
                        }
                        _ => {}
                    }
                }
            }
        }

        marked_nodes.insert(1); // mark root node
        Ok(NCBITaxonomy {
            parent_map,
            child_map,
            rank_map,
            known_ranks,
            name_map,
            marked_nodes,
        })
    }

    pub fn mark_node(&mut self, mut taxid: u64) {
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            taxid = self.parent_map[&taxid];
        }
    }

    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> Result<()> {
        let zero_node = TaxonomyNode::default();
        let node_count = self.marked_nodes.len() as u64 + 1;
        let mut nodes = vec![zero_node.clone(); node_count as usize];

        let mut name_data = Vec::new();
        let mut rank_data = Vec::new();
        let mut rank_offsets = HashMap::new();

        // store ranks
        {
            let mut offset = 0u64;
            for rank in &self.known_ranks {
                rank_offsets.insert(rank.clone(), offset);
                let bytes = rank.as_bytes();
                rank_data.extend_from_slice(bytes);
                rank_data.push(0); // null terminator
                offset += (bytes.len() as u64) + 1;
            }
        }

        let mut external_id_map = HashMap::new();
        external_id_map.insert(0u64, 0u64);
        external_id_map.insert(1u64, 1u64);

        // BFS
        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(1u64);

        let mut internal_node_id = 0u64;
        while let Some(external_node_id) = bfs_queue.pop_front() {
            internal_node_id += 1;
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = zero_node.clone();
            let parent_id = *self.parent_map.get(&external_node_id).unwrap_or(&0);
            node.parent_id = *external_id_map.get(&parent_id).unwrap_or(&0);
            node.external_id = external_node_id;
            node.rank_offset = *rank_offsets.get(&self.rank_map[&external_node_id]).unwrap();
            node.name_offset = name_data.len() as u64;
            // children
            let children = self.child_map.get(&external_node_id);
            let mut child_count = 0;
            // Compute first_child: internal_node_id + 1 + BFS size is a heuristic
            // from original code. But original code sets first_child = internal_node_id + 1 + queue.size()
            // We'll replicate logic: first_child = internal_node_id + 1 + bfs_queue.len()
            node.first_child = internal_node_id + 1 + (bfs_queue.len() as u64);

            if let Some(child_set) = children {
                for &child_node in child_set {
                    if self.marked_nodes.contains(&child_node) {
                        bfs_queue.push_back(child_node);
                        child_count += 1;
                    }
                }
            }
            node.child_count = child_count;
            nodes[internal_node_id as usize] = node;

            let node_name = self
                .name_map
                .get(&external_node_id)
                .cloned()
                .unwrap_or_default();
            name_data.extend_from_slice(node_name.as_bytes());
            name_data.push(0);
        }

        let rank_data_len = rank_data.len() as u64;
        let name_data_len = name_data.len() as u64;

        let taxo = Taxonomy {
            node_count,
            nodes: nodes,
            name_data,
            rank_data,
            name_data_len,
            rank_data_len,
            external_to_internal_id_map: HashMap::new(),
            file_backed: false,
        };

        taxo.write_to_disk(filename)
    }
}

pub struct Taxonomy {
    pub node_count: u64,
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: Vec<u8>,
    pub rank_data: Vec<u8>,
    pub name_data_len: u64,
    pub rank_data_len: u64,
    pub external_to_internal_id_map: HashMap<u64, u64>,
    pub file_backed: bool,
}

impl Taxonomy {
    pub fn new(filename: &str, memory_mapping: bool) -> Result<Taxonomy> {
        // For simplicity, we ignore memory_mapping here
        // and always load into memory.
        // If memory mapping is needed, use memmap2 crate and map the file.
        let mut file = File::open(filename)
            .with_context(|| format!("unable to open taxonomy file {}", filename))?;

        let mut magic_buf = vec![0u8; FILE_MAGIC.len()];
        file.read_exact(&mut magic_buf)?;
        let magic_str = str::from_utf8(&magic_buf).unwrap();
        if magic_str != FILE_MAGIC {
            bail!("attempt to load taxonomy from malformed file {}", filename);
        }

        let mut node_count = 0u64;
        file.read_exact(bytemuck::bytes_of_mut(&mut node_count))?;

        let mut name_data_len = 0u64;
        file.read_exact(bytemuck::bytes_of_mut(&mut name_data_len))?;

        let mut rank_data_len = 0u64;
        file.read_exact(bytemuck::bytes_of_mut(&mut rank_data_len))?;

        let mut nodes = vec![TaxonomyNode::default(); node_count as usize];
        {
            let mut node_bytes =
                vec![0u8; std::mem::size_of::<TaxonomyNode>() * (node_count as usize)];
            file.read_exact(&mut node_bytes)?;
            nodes = bincode::deserialize(&node_bytes)?;
        }

        let mut name_data = vec![0u8; name_data_len as usize];
        file.read_exact(&mut name_data)?;

        let mut rank_data = vec![0u8; rank_data_len as usize];
        file.read_exact(&mut rank_data)?;

        // If not enough data, fail
        // The read_exact calls already ensure completeness.

        Ok(Taxonomy {
            node_count,
            nodes,
            name_data,
            rank_data,
            name_data_len,
            rank_data_len,
            external_to_internal_id_map: HashMap::new(),
            file_backed: false,
        })
    }

    pub fn write_to_disk(&self, filename: &str) -> Result<()> {
        let mut taxo_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(filename)
            .with_context(|| format!("error writing taxonomy to {}", filename))?;

        taxo_file.write_all(FILE_MAGIC.as_bytes())?;
        taxo_file.write_all(bytemuck::bytes_of(&self.node_count))?;
        taxo_file.write_all(bytemuck::bytes_of(&self.name_data_len))?;
        taxo_file.write_all(bytemuck::bytes_of(&self.rank_data_len))?;
        unsafe {
            let node_bytes = std::slice::from_raw_parts(
                self.nodes.as_ptr() as *const u8,
                std::mem::size_of::<TaxonomyNode>() * (self.node_count as usize),
            );
            taxo_file.write_all(node_bytes)?;
        }
        taxo_file.write_all(&self.name_data)?;
        taxo_file.write_all(&self.rank_data)?;

        Ok(())
    }

    pub fn is_a_ancestor_of_b(&self, a: u64, mut b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        while b > a {
            b = self.nodes[b as usize].parent_id;
        }
        b == a
    }

    pub fn lowest_common_ancestor(&self, mut a: u64, mut b: u64) -> u64 {
        if a == 0 || b == 0 {
            // LCA(x,0)=x, LCA(0,x)=x
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

    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.insert(0, 0);
        for i in 1..(self.node_count as usize) {
            let ext_id = self.nodes[i].external_id;
            self.external_to_internal_id_map.insert(ext_id, i as u64);
        }
    }

    pub fn node_count(&self) -> u64 {
        self.node_count
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

    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        *self
            .external_to_internal_id_map
            .get(&external_id)
            .unwrap_or(&0)
    }

    pub fn name_data_at(&self, offset: u64) -> &str {
        let start = offset as usize;
        let end = self.name_data[start..]
            .iter()
            .position(|&c| c == 0)
            .unwrap();
        let slice = &self.name_data[start..start + end];
        std::str::from_utf8(slice).unwrap_or("")
    }

    pub fn rank_data_at(&self, offset: u64) -> &str {
        let start = offset as usize;
        let end = self.rank_data[start..]
            .iter()
            .position(|&c| c == 0)
            .unwrap();
        let slice = &self.rank_data[start..start + end];
        std::str::from_utf8(slice).unwrap_or("")
    }
}
