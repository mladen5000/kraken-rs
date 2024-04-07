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

impl NCBITaxonomy {
    // Existing implementations...

    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(filename)?;
        let mut nodes: Vec<TaxonomyNode> = Vec::new();
        let mut name_data = Vec::new();
        let mut rank_data = Vec::new();

        // Example: Simplified logic for creating TaxonomyNode instances
        for (&external_id, &parent_id) in &self.parent_map {
            if self.marked_nodes.contains(&external_id) {
                let node = TaxonomyNode {
                    parent_id,
                    first_child: 0, // Placeholder, needs to be calculated
                    child_count: 0, // Placeholder, needs to be calculated
                    name_offset: 0, // Placeholder, needs to be calculated
                    rank_offset: 0, // Placeholder, needs to be calculated
                    external_id,
                    godparent_id: 0, // Reserved for future use
                };
                nodes.push(node);
            }
        }

        // Example: Writing the nodes to disk (simplified)
        for node in &nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            // Repeat for other fields...
        }

        Ok(())
    }
}


pub struct Taxonomy {
    file_backed: bool,
    nodes: Vec<TaxonomyNode>,
    name_data: Vec<u8>,
    rank_data: Vec<u8>,
    node_count: usize,
    external_to_internal_id_map: HashMap<u64, u64>,
}

impl Taxonomy {
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