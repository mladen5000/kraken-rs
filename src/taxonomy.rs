use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::error::Error;
use std::fs::File;
use std::hash::Hash;
use std::io::BufRead;
use std::io::BufReader;
use std::io::ErrorKind;
use std::io::Read;
use std::io::Write;
use std::mem;
use std::process;
use std::str::FromStr;

const FILE_MAGIC: &str = "K2TAXDAT";
#[derive(Debug, Default)]
pub struct NCBITaxonomy {
    parent_map_: HashMap<usize, usize>,
    name_map_: HashMap<usize, String>,
    rank_map_: HashMap<usize, String>,
    child_map_: HashMap<usize, HashSet<usize>>,
    marked_nodes_: HashSet<usize>,
    known_ranks_: HashSet<String>,
}

#[derive(Clone, Debu)]
struct Taxonomy {
    nodes_: Vec<TaxonomyNode>,
    node_count_: usize,
    name_data_: String,
    name_data_len_: usize,
    rank_data_: String,
    rank_data_len_: usize,
    external_to_internal_id_map_: HashMap<usize, usize>,
}
impl Default for Taxonomy {
    fn default() -> Self {
        Taxonomy {
            ..Default::default()
        }
    }
}

impl Taxonomy {
    fn init(&mut self, filename: &str, memory_mapping: bool) -> Result<(), std::io::Error> {
        if memory_mapping {
            // Memory mapping is not directly supported in Rust's standard library.
            // You would need to use a crate like memmap.
            unimplemented!()
        } else {
            let mut file = File::open(filename)?;

            let mut magic = vec![0; FILE_MAGIC.len()];
            file.read_exact(&mut magic)?;
            if magic != FILE_MAGIC.as_bytes() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "malformed taxonomy file",
                ));
            }

            let mut buffer = [0; mem::size_of::<usize>()];
            file.read_exact(&mut buffer)?;
            self.node_count_ = usize::from_le_bytes(buffer);

            file.read_exact(&mut buffer)?;
            self.name_data_len_ = usize::from_le_bytes(buffer);

            file.read_exact(&mut buffer)?;
            self.rank_data_len_ = usize::from_le_bytes(buffer);

            self.nodes_ = vec![TaxonomyNode::default(); self.node_count_ as usize];
            for node in &mut self.nodes_ {
                file.read_exact(&mut buffer)?;
                node.parent_id = usize::from_le_bytes(buffer);
            }

            let mut name_data = vec![0; self.name_data_len_ as usize];
            file.read_exact(&mut name_data)?;
            self.name_data_ = String::from_utf8_lossy(&name_data).to_string();

            let mut rank_data_ = vec![0; self.rank_data_len_ as usize];
            file.read_exact(&mut rank_data_)?;
            self.rank_data_ = String::from_utf8_lossy(&rank_data_).to_string();
        }
        Ok(())
    }
    /// Logic here depends on higher nodes having smaller IDs
    /// Idea: advance B tracker up tree, A is ancestor iff B tracker hits A
    fn is_a_ancestor_of_b(&self, mut a: usize, mut b: usize) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        while b > a {
            b = self.nodes_[b].parent_id;
        }
        b == a
    }
    /// Logic here depends on higher nodes having smaller IDs
    /// Idea: track two nodes, advance lower tracker up tree, trackers meet @ LCA
    fn lowest_common_ancestor(&self, mut a: usize, mut b: usize) -> usize {
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }
        while a != b {
            if a > b {
                a = self.nodes_[a].parent_id;
            } else {
                b = self.nodes_[b].parent_id;
            }
        }
        a
    }

    fn write_to_disk(&self, filename: &str) {
        let mut taxo_file = File::create(filename).expect("Unable to create file");
        taxo_file
            .write_all(FILE_MAGIC.as_bytes())
            .expect("Unable to write data");
        taxo_file
            .write_all(&self.node_count_.to_le_bytes())
            .expect("Unable to write data");
        taxo_file
            .write_all(&self.name_data_len_.to_le_bytes())
            .expect("Unable to write data");
        taxo_file
            .write_all(&self.rank_data_len_.to_le_bytes())
            .expect("Unable to write data");
        for node in &self.nodes_ {
            taxo_file
                .write_all(unsafe {
                    std::slice::from_raw_parts(
                        (node as *const _ as *const u8),
                        mem::size_of_val(node),
                    )
                })
                .expect("Unable to write data");
        }

        taxo_file
            .write_all(self.name_data_.as_bytes())
            .expect("Unable to write data");
        taxo_file
            .write_all(self.rank_data_.as_bytes())
            .expect("Unable to write data");
    }
}

struct TaxonomyNode {
    parent_id: usize,
    first_child: usize,
    child_count: usize,
    name_offset: usize,
    rank_offset: usize,
    external_id: usize,
    godparent_id: usize,
}
impl Default for TaxonomyNode {
    fn default() -> Self {
        TaxonomyNode {
            first_child: 0,
            parent_id: 0,
            child_count: 0,
            name_offset: 0,
            rank_offset: 0,
            external_id: 0,
            godparent_id: 0,
        }
    }
}
impl Clone for TaxonomyNode {
    fn clone(&self) -> Self {
        TaxonomyNode {
            first_child: self.first_child,
            parent_id: self.parent_id,
            child_count: self.child_count,
            name_offset: self.name_offset,
            rank_offset: self.name_offset,
            external_id: self.external_id,
            godparent_id: self.godparent_id,
        }
    }
}

impl NCBITaxonomy {
    pub fn new(&self, nodes_filename: &str, names_filename: &str) -> Self {
        // Open nodes and names
        let nodes_file = get_bufreader(nodes_filename);
        let names_file = get_bufreader(names_filename);

        // todo: delete these
        let line: String;
        let (mut node_id, mut parent_id): (usize, usize);
        let (mut name, mut rank, junk): (String, String, String);
        const DELIM: &str = "\t|\t";

        while let Some(line) = nodes_file.lines().next() {
            let mut line = line.unwrap();
            line.pop(); // discard trailing
            line.pop(); //   "\t|"
            let mut pos1 = 0usize;
            let mut field_ct: i64 = 0;
            let mut finished = false;
            let mut token;
            while field_ct < 10 && !finished {
                field_ct += 1;
                let pos2 = line.find(DELIM);
                if pos2 == None {
                    token = line.get(pos1..);
                    finished = true
                } else {
                    token = line.get(pos1..(pos2.unwrap() - pos1));
                    pos1 = pos2.unwrap() + DELIM.len();
                }

                match field_ct {
                    // 1-based counting
                    1 => {
                        node_id = token.unwrap().parse::<usize>().unwrap();
                        if node_id == 0 {
                            panic!("Attempt to create taxonomy with node ID == 0")
                        }
                    }
                    2 => parent_id = token.unwrap().parse::<usize>().unwrap(),
                    3 => {
                        rank = token.unwrap().to_owned();
                        finished = true;
                    }
                    _ => {}
                }
                if node_id == 1 {
                    parent_id = 0;
                }

                self.parent_map_.insert(node_id, parent_id);
                self.child_map_
                    .entry(parent_id)
                    .or_insert_with(HashSet::new)
                    .insert(node_id);
                self.rank_map_.insert(node_id, rank);
                self.known_ranks_.insert(rank);
            }
        }

        while let Some(line) = names_file.lines().next() {
            let mut line = line.unwrap();
            line.pop(); // discard trailing
            line.pop(); //   "\t|"
            let mut pos1 = 0usize;
            let mut field_ct = 0i64;
            let mut finished = false;
            let mut token;
            ////
            while field_ct < 10 && !finished {
                let pos2 = line.find(DELIM);
                // let token: String;
                if pos2 == None {
                    token = line.get(pos1..);
                    finished = true
                } else {
                    token = line.get(pos1..(pos2.unwrap() - pos1));
                    pos1 = pos2.unwrap() + DELIM.len();
                }
                match field_ct {
                    // 1-based counting
                    // ...
                    1 => {
                        node_id = token.unwrap().parse::<usize>().unwrap();
                        if node_id == 0 {
                            panic!("Attempt to create taxonomy with node ID == 0")
                        }
                    }
                    2 => name = token.unwrap().to_owned(),
                    4 => {
                        if token == Some("scientific name") {
                            self.name_map_.insert(node_id, name);
                        }
                        finished = true
                    }
                    _ => {}
                }
            }
        }
        self.marked_nodes_.insert(1);
        *self
    }

    pub fn mark_node(&self, mut taxid: usize) {
        while !self.marked_nodes_.contains(&taxid) {
            self.marked_nodes_.insert(taxid);
            taxid = self.parent_map_[&taxid]
        }
    }

    fn convert_to_kraken_taxonomy(&self, filename: &str) {
        let mut taxo = Taxonomy::default();
        let zeroes_node = TaxonomyNode::default();

        taxo.node_count_ = self.marked_nodes_.len() + 1; // +1 because 0 is illegal
        taxo.nodes_ = vec![TaxonomyNode::default(); taxo.node_count_];
        let mut name_data: String;
        let mut rank_data: String;

        // Because so many of the node rank names are shared, we only store one copy of each rank
        let mut rank_offsets: HashMap<String, usize> = HashMap::new();
        for rank in self.known_ranks_ {
            rank_offsets.insert(rank, rank_data.len());
            rank_data.push_str(&rank);
            rank_data.push('\0');
        }
        let mut internal_node_id = 0;
        let mut external_id_map: HashMap<usize, usize> = HashMap::new();
        external_id_map[&0] = 0;
        external_id_map[&1] = 1; // 1 is root in both NCBI and Kraken taxonomies

        // Breadth-first search (BFS) through NCBI taxonomy, assigning internal IDs
        // In sequential order as nodes are encountered via BFS.
        let mut bfs_queue: VecDeque<usize>;
        bfs_queue.push_back(1);
        while !bfs_queue.is_empty() {
            internal_node_id += 1;
            let external_node_id = bfs_queue.front().unwrap(); // Unwrap the Option<&usize> value
            bfs_queue.pop_back();
            external_id_map.insert(*external_node_id, internal_node_id); // Fix: Dereference external_node_id

            //todo: Probably could just be initialized
            let mut node = zeroes_node;
            node.parent_id = external_id_map[&self.parent_map_[external_node_id]]; // Fix: Dereference external_node_id
            node.external_id = *external_node_id; // Fix: Dereference external_node_id
            node.rank_offset = rank_offsets[&self.rank_map_[external_node_id]]; // Fix: Dereference external_node_id
            node.name_offset = name_data.len();
            node.first_child = internal_node_id + 1 + bfs_queue.len();
            for child_node in self.child_map_[external_node_id] {
                if self.marked_nodes_.contains(&child_node) {
                    bfs_queue.push_back(child_node);
                    node.child_count += 1
                }
            }
            taxo.nodes_[internal_node_id] = node;
            let name: String = self.name_map_[external_node_id];
            name_data.push_str(&name)
        } // end BFS while loop

        taxo.rank_data_ = rank_data;
        taxo.rank_data_len_ = rank_data.len();

        taxo.name_data_ = name_data;
        taxo.name_data_len_ = name_data.len();

        taxo.write_to_disk(filename);
    }
}

fn get_bufreader(filename: &str) -> BufReader<File> {
    match File::open(filename) {
        Ok(file) => BufReader::new(file),
        Err(_) => {
            eprintln!("Cannot open {}", filename);
            process::exit(1);
        }
    }
}