#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ncbi_taxonomy_new() {
        let ncbi_taxonomy = NCBITaxonomy::new("path/to/nodes.dmp", "path/to/names.dmp").unwrap();
        assert!(!ncbi_taxonomy.parent_map.is_empty());
        assert!(!ncbi_taxonomy.name_map.is_empty());
        assert!(!ncbi_taxonomy.rank_map.is_empty());
        assert!(!ncbi_taxonomy.child_map.is_empty());
        assert!(ncbi_taxonomy.marked_nodes.contains(&1)); // Root node is marked
    }

    #[test]
    fn test_convert_to_kraken_taxonomy() {
        let ncbi_taxonomy = NCBITaxonomy::new("path/to/nodes.dmp", "path/to/names.dmp").unwrap();
        let result = ncbi_taxonomy.convert_to_kraken_taxonomy("path/to/output.k2t");
        assert!(result.is_ok());
        // Additional assertions to verify the contents of the output file...
    }
