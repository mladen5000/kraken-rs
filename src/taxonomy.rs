use std::collections::{HashMap, HashSet};

#[derive(Debug)]
pub struct TaxonomyNode {
    parent_id: u64,
    first_child: u64,
    child_count: u64,
    name_offset: u64,
    rank_offset: u64,
    external_id: u64,
    godparent_id: u64,
}

pub struct NCBITaxonomy {
    nodes_filename: String,
    names_filename: String,
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::error::Error;

    pub fn new(nodes_filename: &str, names_filename: &str) -> Result<Self, Box<dyn Error>> {
        let nodes_file = File::open(nodes_filename)
            .map_err(|e| format!("Failed to open nodes file '{}': {}", nodes_filename, e))?;
        let names_file = File::open(names_filename)
            .map_err(|e| format!("Failed to open names file '{}': {}", names_filename, e))?;

        // Continue with the rest of the method...
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
}