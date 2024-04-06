/*
* Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
*
* This file is part of the Kraken 2 taxonomic sequence classification system.
*/

use std::fs::File;
use std::io;

use memmap2::MmapMut;

pub struct MMapFile {
    mmap: Option<MmapMut>,
    file: Option<File>,
    filesize: usize,
}

impl MMapFile {
    pub fn new() -> Self {
        MMapFile {
            mmap: None,
            file: None,
            filesize: 0,
        }
    }

    pub fn open_file(&mut self, filename: &str, size: usize) -> io::Result<()> {
        let file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .open(filename)?;

        if size > 0 {
            file.set_len(size as u64)?;
            self.filesize = size;
        } else {
            self.filesize = file.metadata()?.len() as usize;
        }

        let mmap = unsafe { MmapMut::map_mut(&file)? };
        self.file = Some(file);
        self.mmap = Some(mmap);
        Ok(())
    }

    pub fn load_file(&mut self) {
        if let Some(mmap) = &self.mmap {
            let page_size = 4096;

            for pos in (0..self.filesize).step_by(page_size) {
                let end = std::cmp::min(pos + page_size, self.filesize);
                let len = end - pos;
                let mut buf = vec![0u8; len];
                buf.copy_from_slice(&mmap.as_ref()[pos..end]);
            }
        }
    }

    pub fn fptr(&self) -> *const u8 {
        match &self.mmap {
            Some(mmap) => mmap.as_ptr(),
            None => std::ptr::null(),
        }
    }

    pub fn filesize(&self) -> usize {
        self.filesize
    }

    pub fn sync_file(&mut self) -> io::Result<()> {
        if let Some(mmap) = &mut self.mmap {
            mmap.flush()?;
        }
        Ok(())
    }

    pub fn close_file(&mut self) -> io::Result<()> {
        self.sync_file()?;
        self.mmap = None;
        if let Some(file) = &self.file {
            file.sync_all()?;
        }
        self.file = None;
        Ok(())
    }
}

impl Drop for MMapFile {
    fn drop(&mut self) {
        let _ = self.close_file();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[test]
    fn test_open_and_close_file() {
        let filename = "test_open_and_close_file.txt";
        let mut file = File::create(filename).unwrap();
        file.write_all(b"Hello, world!").unwrap();

        let mut mmap_file = MMapFile::new();
        mmap_file.open_file(filename, 0).unwrap();
        assert!(mmap_file.mmap.is_some());
        assert!(mmap_file.file.is_some());
        assert!(mmap_file.filesize() > 0);

        mmap_file.close_file().unwrap();
        assert!(mmap_file.mmap.is_none());
        assert!(mmap_file.file.is_none());

        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn test_load_file() {
        let filename = "test_load_file.txt";
        let mut file = File::create(filename).unwrap();
        let data = "Hello, world!";
        file.write_all(data.as_bytes()).unwrap();

        let mut mmap_file = MMapFile::new();
        mmap_file.open_file(filename, 0).unwrap();
        mmap_file.load_file();

        let loaded_data = unsafe { std::slice::from_raw_parts(mmap_file.fptr(), data.len()) };
        assert_eq!(loaded_data, data.as_bytes());

        mmap_file.close_file().unwrap();
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn test_sync_file() {
        let filename = "test_sync_file.txt";
        let mut file = File::create(filename).unwrap();
        let data = "Hello, world!";
        file.write_all(data.as_bytes()).unwrap();

        let mut mmap_file = MMapFile::new();
        mmap_file.open_file(filename, 0).unwrap();
        let fptr = mmap_file.fptr();
        unsafe {
            let slice = std::slice::from_raw_parts_mut(fptr as *mut u8, data.len());
            slice.copy_from_slice(b"Modified data");
        }
        mmap_file.sync_file().unwrap();
        mmap_file.close_file().unwrap();

        let modified_data = fs::read_to_string(filename).unwrap();
        assert_eq!(modified_data, "Modified data");

        fs::remove_file(filename).unwrap();
    }
}
