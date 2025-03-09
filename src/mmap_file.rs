/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 * Rust conversion Copyright 2025
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use anyhow::{Context, Result};
use memmap2::MmapOptions;
use std::fs::File;
use std::path::Path;

/// A memory-mapped file, equivalent to the MMapFile class in C++.
pub struct MMapFile {
    file: File,
    mmap: memmap2::Mmap,
}

impl MMapFile {
    /// Opens a file for memory mapping.
    pub fn open_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("Failed to open file for memory mapping: {:?}", path.as_ref()))?;
        
        let mmap = unsafe { 
            MmapOptions::new()
                .map(&file)
                .with_context(|| format!("Failed to memory map file: {:?}", path.as_ref()))?
        };
        
        Ok(MMapFile { file, mmap })
    }

    /// Returns a pointer to the start of the memory-mapped data.
    pub fn fptr(&self) -> &[u8] {
        &self.mmap
    }

    /// Returns the size of the memory-mapped file in bytes.
    pub fn len(&self) -> usize {
        self.mmap.len()
    }
    
    /// Returns true if the memory-mapped file is empty.
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
    
    /// Returns the file size.
    pub fn filesize(&self) -> usize {
        self.mmap.len()
    }
}