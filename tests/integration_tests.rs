#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn test_taxonomy_conversion() -> Result<(), Box<dyn std::error::Error>> {
        let nodes_filename = "tests/data/nodes.dmp";
        let names_filename = "tests/data/names.dmp";
        let output_filename = "tests/output/output.k2t";

        // Ensure the output directory exists
        fs::create_dir_all("tests/output")?;

        let ncbi_taxonomy = NCBITaxonomy::new(nodes_filename, names_filename)?;
        ncbi_taxonomy.convert_to_kraken_taxonomy(output_filename)?;

        // Verify the output file exists
        assert!(Path::new(output_filename).exists());

        // Additional verification steps could include checking the contents of the output file
        // to ensure it matches the expected format and data for the Kraken taxonomy.

        Ok(())
    }
}

