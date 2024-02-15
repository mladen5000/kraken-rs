use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use std::mem;
use std::process;

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
impl NCBITaxonomy {
    pub fn new(&mut self, nodes_filename: &str, names_filename: &str) -> Self {
        const DELIM: &str = "\t|\t";
        // Open nodes and names
        let nodes_file = get_bufreader(nodes_filename);
        let names_file = get_bufreader(names_filename);

        // todo: delete these
        let line: String;
        let (mut node_id, mut parent_id): (usize, usize) = (0, 0);
        let mut rank = String::new();
        let mut name = String::new();
        let mut junk = String::new();

        let mut lines = nodes_file.lines();
        while let Some(line) = lines.next() {
            let mut line = line.unwrap();
            line.pop(); // discard trailing
            line.pop(); //   "\t|"
            let mut pos1 = 0usize;
            let mut field_ct = 0i64;
            let mut finished = false;
            let mut token;
            while field_ct < 10 && !finished {
                field_ct += 1;
                match line.find(DELIM) {
                    Some(pos2) => {
                        token = line.get(pos1..(pos2 - pos1));
                        pos1 = pos2 + DELIM.len();
                    }
                    None => {
                        token = line.get(pos1..);
                        finished = true;
                    }
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
                        let rank = token.unwrap().to_owned();
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
                self.rank_map_.insert(node_id, rank.to_owned());
                self.known_ranks_.insert(rank.to_string());
            }
        }

        let mut lines = names_file.lines();
        while let Some(line) = lines.next() {
            let mut line = line.unwrap();

            line.pop(); // discard trailing
            line.pop(); //   "\t|"

            let mut pos1 = 0usize;
            let mut field_ct = 0;
            let mut finished = false;
            let mut token: Option<&str> = None;
            while field_ct < 10 && !finished {
                field_ct += 1;
                let name = name.clone();
                // Rest of the code...
            }
            match line.find(DELIM) {
                Some(pos2) => {
                    token = Some(&line[pos1..pos2]);
                    pos1 = pos2 + DELIM.len();
                }
                None => {
                    token = Some(&line);
                    finished = true;
                }
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
                2 => {
                    name = token.unwrap().to_owned();
                }
                4 => {
                    if token == Some("scientific name") {
                        self.name_map_.insert(node_id, name);
                    }
                    finished = true
                }
                // Remove the closing curly brace and the empty block
                _ => todo!(),
            }
        }
        self.marked_nodes_.insert(1);
        *self
    }
    pub fn mark_node(mut self, mut taxid: usize) -> Self {
        while !self.marked_nodes_.contains(&taxid) {
            self.marked_nodes_.insert(taxid);
            taxid = self.parent_map_[&taxid]
        }
        self
    }
    pub fn convert_to_kraken_taxonomy(&self, filename: &str) {
        let mut taxo = Taxonomy::default();
        let zeroes_node = TaxonomyNode::default().clone();

        taxo.node_count_ = self.marked_nodes_.len() + 1; // +1 because 0 is illegal
        taxo.nodes_ = vec![TaxonomyNode::default(); taxo.node_count_];
        let mut name_data: String;
        let mut rank_data: String;

        // Because so many of the node rank names are shared, we only store one copy of each rank
        let mut rank_offsets: HashMap<String, usize> = HashMap::new();
        let mut internal_node_id = 0;
        let mut external_id_map: HashMap<usize, usize> = HashMap::new();
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1); // 1 is root in both NCBI and Kraken taxonomies

        // Breadth-first search (BFS) through NCBI taxonomy, assigning internal IDs
        // In sequential order as nodes are encountered via BFS.
        let mut bfs_queue: VecDeque<usize> = VecDeque::new();
        bfs_queue.push_back(1);
        let mut node = zeroes_node; // Move the initialization of `node` outside of the loop
        while !bfs_queue.is_empty() {
            let mut node = zeroes_node.clone();
            internal_node_id += 1;
            let external_node_id = bfs_queue.front().unwrap().clone(); // Clone the Option<&usize> value
            bfs_queue.pop_back();
            external_id_map.insert(external_node_id, internal_node_id); // Remove the dereference operator (*)

            //todo: Probably could just be initialized
            let mut node = zeroes_node;
            node.parent_id = external_id_map[&self.parent_map_[&external_node_id]];
            node.external_id = external_node_id;
            node.rank_offset = rank_offsets[&self.rank_map_[&external_node_id]];
            node.name_offset = name_data.len();
            node.first_child = internal_node_id + 1 + bfs_queue.len();
            for child_node in &self.child_map_[&external_node_id] {
                if self.marked_nodes_.contains(child_node) {
                    bfs_queue.push_back(*child_node);
                    node.child_count += 1;
                }
            }
            taxo.nodes_[internal_node_id] = node;
            let name: &String = &self.name_map_[&external_node_id];
            name_data.push_str(&name)
        } // end BFS while loop

        taxo.rank_data_ = rank_data;
        taxo.rank_data_len_ = rank_data.len();

        taxo.name_data_ = name_data;
        taxo.name_data_len_ = name_data.len();

        taxo.write_to_disk(filename);
    }
}

#[derive(Clone, Default)]
pub struct Taxonomy {
    pub nodes_: Vec<TaxonomyNode>,
    pub node_count_: usize,
    pub name_data_: String,
    pub name_data_len_: usize,
    pub rank_data_: String,
    pub rank_data_len_: usize,
    pub external_to_internal_id_map_: HashMap<usize, usize>,
}

impl Taxonomy {
    pub fn new(filename: &str, memory_mapping: bool) -> Result<Self, std::io::Error> {
        let mut taxonomy = Taxonomy {
            nodes_: Vec::new(),
            node_count_: 0,
            name_data_: String::new(),
            name_data_len_: 0,
            rank_data_: String::new(),
            rank_data_len_: 0,
            external_to_internal_id_map_: HashMap::new(),
        };

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
            taxonomy.node_count_ = usize::from_le_bytes(buffer);

            file.read_exact(&mut buffer)?;
            taxonomy.name_data_len_ = usize::from_le_bytes(buffer);

            file.read_exact(&mut buffer)?;
            taxonomy.rank_data_len_ = usize::from_le_bytes(buffer);

            taxonomy.nodes_ = vec![TaxonomyNode::default(); taxonomy.node_count_ as usize];
            for node in &mut taxonomy.nodes_ {
                file.read_exact(&mut buffer)?;
                node.parent_id = usize::from_le_bytes(buffer);
            }

            let mut name_data = vec![0; taxonomy.name_data_len_ as usize];
            file.read_exact(&mut name_data)?;
            taxonomy.name_data_ = String::from_utf8_lossy(&name_data).to_string();

            let mut rank_data_ = vec![0; taxonomy.rank_data_len_ as usize];
            file.read_exact(&mut rank_data_)?;
            taxonomy.rank_data_ = String::from_utf8_lossy(&rank_data_).to_string();
        }

        Ok(taxonomy)
    }
    pub fn name_data(&self) -> &String {
        &self.name_data_
    }
    pub fn rank_data(&self) -> &String {
        &self.rank_data_
    }
    pub fn node_count(&self) -> usize {
        self.node_count_
    }
    pub fn nodes(&self) -> &Vec<TaxonomyNode> {
        &self.nodes_
    }
    pub fn is_a_ancestor_of_b(&self, mut a: usize, mut b: usize) -> bool {
        // Idea: advance B tracker up tree, A is ancestor iff B tracker hits A
        if a == 0 || b == 0 {
            return false;
        }
        while b > a {
            b = self.nodes_[b].parent_id;
        }
        b == a
    }
    pub fn lowest_common_ancestor(&self, mut a: usize, mut b: usize) -> usize {
        // Logic here depends on higher nodes having smaller IDs
        // Idea: track two nodes, advance lower tracker up tree, trackers meet @ LCA
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
    pub fn write_to_disk(&self, filename: &str) {
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

pub struct TaxonomyNode {
    /// Must be lower-numbered node
    pub parent_id: u64,
    /// Must be higher-numbered node
    pub first_child: u64,
    /// Children of a node are in contiguous block
    pub child_count: u64,
    /// Location of name in name data super-string
    pub name_offset: u64,
    /// Location of rank in rank data super-string
    pub rank_offset: u64,
    /// Taxonomy ID for reporting purposes (usually NCBI)
    pub external_id: u64,
    /// Reserved for future use to enable faster traversal
    pub godparent_id: u64,
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

pub fn get_bufreader(filename: &str) -> BufReader<File> {
    match File::open(filename) {
        Ok(file) => BufReader::new(file),
        Err(_) => {
            eprintln!("Cannot open {}", filename);
            process::exit(1);
        }
    }
}
