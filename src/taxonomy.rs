use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug)]
struct TaxonomyNode {
    parent_id: u64,
    first_child: u64,
    child_count: u64,
    name_offset: u64,
    rank_offset: u64,
    external_id: u64,
    godparent_id: u64,
}

struct NCBITaxonomy {
    nodes_filename: String,
    names_filename: String,
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

impl Taxonomy {
    // Existing implementations...

    pub fn write_to_disk(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(filename)?;

        // Write the file magic header
        file.write_all(b"K2TAXDAT")?;

        // Serialize and write node_count, name_data_len, and rank_data_len
        file.write_all(&self.node_count.to_le_bytes())?;
        file.write_all(&(self.name_data.len() as u64).to_le_bytes())?;
        file.write_all(&(self.rank_data.len() as u64).to_le_bytes())?;

        // Serialize and write each TaxonomyNode
        for node in &self.nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            file.write_all(&node.first_child.to_le_bytes())?;
            file.write_all(&node.child_count.to_le_bytes())?;
            file.write_all(&node.name_offset.to_le_bytes())?;
            file.write_all(&node.rank_offset.to_le_bytes())?;
            file.write_all(&node.external_id.to_le_bytes())?;
            file.write_all(&node.godparent_id.to_le_bytes())?;
        }

        // Write name_data and rank_data
        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        Ok(())
    }
}