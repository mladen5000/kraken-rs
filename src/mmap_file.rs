use std::fs::{File, OpenOptions};
use std::io::{self, Seek, SeekFrom};
use std::os::unix::io::AsRawFd;
use std::os::unix::prelude::FileExt;
use std::path::Path;
use std::sync::Arc;

use memmap2::MmapMut;

pub struct MMapFile {
    mmap: Option<MmapMut>,
    file: Option<File>,
}

impl MMapFile {
    pub fn new() -> Self {
        Self {
            mmap: None,
            file: None,
        }
    }

    pub fn open_file<P: AsRef<Path>>(&mut self, path: P, size: usize) -> io::Result<()> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(path)?;

        file.set_len(size as u64)?;

        let mmap = unsafe { MmapMut::map_mut(&file)? };

        self.file = Some(file);
        self.mmap = Some(mmap);

        Ok(())
    }

    pub fn load_file(&mut self) {
        if let Some(mmap) = &self.mmap {
            let data = Arc::new(mmap.clone());
            let mut handles = vec![];

            for chunk in data.chunks(4096) {
                let data = Arc::clone(&data);
                let handle = std::thread::spawn(move || {
                    let mut buf = [0; 4096];
                    buf.copy_from_slice(chunk);
                });
                handles.push(handle);
            }

            for handle in handles {
                handle.join().unwrap();
            }
        }
    }

    pub fn sync_file(&mut self) {
        if let Some(mmap) = &mut self.mmap {
            mmap.flush().expect("Failed to sync mmap file");
        }
    }

    pub fn close_file(&mut self) {
        self.sync_file();
        self.mmap = None;
        self.file = None;
    }
}

impl Drop for MMapFile {
    fn drop(&mut self) {
        self.close_file();
    }
}
