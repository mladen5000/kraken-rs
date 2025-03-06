use memmap2::MmapOptions;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

pub struct MMapFile {
    file: File,
    mmap: memmap2::Mmap,
}

impl MMapFile {
    pub fn open_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(MMapFile { file, mmap })
    }

    pub fn fptr(&self) -> &[u8] {
        &self.mmap
    }

    pub fn len(&self) -> usize {
        self.mmap.len()
    }
}
