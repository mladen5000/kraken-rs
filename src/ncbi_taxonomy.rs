struct TaxonomyNode {
    /// Must be lower-numbered node
    parent_id: u64,
    /// Must be higher-numbered node
    first_child: u64,
    /// Children of a node are in contiguous block
    child_count: u64,
    /// Location of name in name data super-string
    name_offset: u64,
    /// Location of rank in rank data super-string
    rank_offset: u64,
    /// Taxonomy ID for reporting purposes (usually NCBI)
    external_id: u64,
    /// Reserved for future use to enable faster traversal
    godparent_id: u64,
}
