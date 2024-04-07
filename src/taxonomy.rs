// Placeholder for output equivalence testing...

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::error::Error;

pub struct NCBITaxonomy {
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

impl NCBITaxonomy {
    pub fn new(nodes_filename: &str, names_filename: &str) -> Result<Self, Box<dyn Error>> {
        let mut ncbi_taxonomy = NCBITaxonomy {
            parent_map: HashMap::new(),
            name_map: HashMap::new(),
            rank_map: HashMap::new(),
            child_map: HashMap::new(),
            marked_nodes: HashSet::new(),
            known_ranks: HashSet::new(),
        };

        let mut parent_map: HashMap<u64, u64> = HashMap::new();
        let mut name_map: HashMap<u64, String> = HashMap::new();
        let mut rank_map: HashMap<u64, String> = HashMap::new();
        let mut child_map: HashMap<u64, HashSet<u64>> = HashMap::new();
        let mut marked_nodes: HashSet<u64> = HashSet::new();
        let mut known_ranks: HashSet<String> = HashSet::new();

        let nodes_file = File::open(nodes_filename)?;
        let reader = BufReader::new(nodes_file);
        // Parse nodes file
        let nodes_file = File::open(nodes_filename)?;
        let reader = BufReader::new(nodes_file);
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split("\t|\t").collect();
            if parts.len() < 3 { continue; }

            let node_id: u64 = parts[0].parse()?;
            let parent_id: u64 = parts[1].parse()?;
            let rank: String = parts[2].to_string();

            ncbi_taxonomy.parent_map.insert(node_id, parent_id);
            ncbi_taxonomy.rank_map.insert(node_id, rank);
            ncbi_taxonomy.child_map.entry(parent_id).or_insert_with(HashSet::new).insert(node_id);
        }

        // Parse names file
        let names_file = File::open(names_filename)?;
        let names_reader = BufReader::new(names_file);
        for line in names_reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split("\t|\t").collect();
            if parts.len() < 4 { continue; }

            let node_id: u64 = parts[0].parse()?;
            let name: String = parts[1].to_string();

            ncbi_taxonomy.name_map.insert(node_id, name);
        }


        Ok(NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            marked_nodes,
            known_ranks,
        })
    }
}
}
}




    pub fn move_to_memory(&mut self) -> Result<(), Box<dyn Error>> {
        // Assuming the initial loading mechanism involves memory-mapped files or similar,
        // this method would convert those structures into purely in-memory equivalents.

        // Example: Convert memory-mapped name and rank data into Vec<u8>
        self.name_data = Vec::from(self.name_data.as_slice());
        self.rank_data = Vec::from(self.rank_data.as_slice());

        // Convert nodes from a memory-mapped array to a Vec<TaxonomyNode>
        self.nodes = self.nodes.clone();

        Ok(())
    }


    // Existing implementations...

    pub fn is_a_ancestor_of_b(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        let mut current = b;
        while current > a && current != 0 {
            current = self.nodes[current as usize].parent_id;
        }
        current == a
    }

    pub fn lowest_common_ancestor(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            return if a == 0 { b } else { a };
        }
        let mut ancestors_a = HashSet::new();
        let mut current_a = a;
        while current_a != 0 {
            ancestors_a.insert(current_a);
            current_a = self.nodes[current_a as usize].parent_id;
        }
        let mut current_b = b;
        while !ancestors_a.contains(&current_b) && current_b != 0 {
            current_b = self.nodes[current_b as usize].parent_id;
        }
        current_b
    }
    pub fn mark_node(&mut self, taxid: u64) {
        let mut current_taxid = taxid;
        while !self.marked_nodes.contains(&current_taxid) {
            self.marked_nodes.insert(current_taxid);
            if let Some(&parent_id) = self.parent_map.get(&current_taxid) {
                current_taxid = parent_id;
            } else {
                break; // Reached the root or an unconnected node
            }
        }
    }
    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        // Write the header or initial taxonomy data to the file
        writer.write_all(b"Kraken taxonomy data\n")?;

        // Iterate over marked nodes and organize them for Kraken format
        for (&node_id, _) in self.marked_nodes.iter() {
            let parent_id = self.parent_map.get(&node_id).unwrap_or(&0);
            let rank = self.rank_map.get(&node_id).unwrap_or(&"unknown".to_string());
            let name = self.name_map.get(&node_id).unwrap_or(&"unknown".to_string());

            // Writing node data to the file
            writer.write_fmt(format_args!("node_id: {}, parent_id: {}, rank: {}, name: {}\n", node_id, parent_id, rank, name))?;
        }

        writer.flush()?;
        Ok(())
    }
    }
}