use std::fs::{File, OpenOptions};
use std::io;
use std::io::Read;
use std::os::unix::fs::FileExt;

use std::path::Path;
use std::sync::Arc;

use memmap2::Mmap;

pub struct MMapFile {
    mmap: Option<Mmap>,
    file: Option<File>,
}

impl MMapFile {
    pub fn new() -> Self {
        Self {
            mmap: None,
            file: None,
        }
    }
    pub fn fptr(&self) -> *mut u8 {
        match &self.mmap {
            Some(mmap) => mmap.as_ptr() as *mut u8,
            None => std::ptr::null_mut(),
        }
    }

    pub fn filesize(&self) -> usize {
        match &self.mmap {
            Some(mmap) => mmap.len(),
            None => 0,
        }
    }

    pub fn open_file<P: AsRef<Path>>(&mut self, path: P, size: usize) -> io::Result<()> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(path)?;

        file.set_len(size as u64)?;

        let mmap = unsafe { Mmap::map(&file)? };

        self.file = Some(file);
        self.mmap = Some(mmap);

        Ok(())
    }

    pub fn load_file(&self) {
        if let Some(mmap) = &self.mmap {
            let data = mmap.as_ref();
            let mut handles = vec![];

            for chunk in data.chunks(4096) {
                let chunk = chunk.to_vec();
                let handle = std::thread::spawn(move || {
                    let mut buf = [0; 4096];
                    buf[..chunk.len()].copy_from_slice(&chunk);
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
            let ptr = mmap.as_ptr();
            let len = mmap.len();

            unsafe {
                let result = libc::msync(ptr as *mut libc::c_void, len, libc::MS_SYNC);
                if result != 0 {
                    panic!(
                        "Failed to sync mmap file: {}",
                        std::io::Error::last_os_error()
                    );
                }
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::remove_file;
    use std::io::Write;

    #[test]
    fn test_open_and_close_file() {
        let mut mmap_file = MMapFile::new();
        assert!(mmap_file.open_file("test.txt", 100).is_ok());
        assert!(mmap_file.file.is_some());
        assert!(mmap_file.mmap.is_some());
        mmap_file.close_file();
        assert!(mmap_file.file.is_none());
        assert!(mmap_file.mmap.is_none());
        remove_file("test.txt").unwrap();
    }

    #[test]
    fn test_filesize() {
        let mut mmap_file = MMapFile::new();
        mmap_file.open_file("test.txt", 100).unwrap();
        assert_eq!(mmap_file.filesize(), 100);
        mmap_file.close_file();
        remove_file("test.txt").unwrap();
    }

    #[test]
    fn test_fptr() {
        let mut mmap_file = MMapFile::new();
        mmap_file.open_file("test.txt", 100).unwrap();
        let fptr = mmap_file.fptr();
        assert!(!fptr.is_null());
        mmap_file.close_file();
        remove_file("test.txt").unwrap();
    }

    #[test]
    fn test_load_file() {
        let mut mmap_file = MMapFile::new();
        mmap_file.open_file("test.txt", 100).unwrap();
        let data = b"Hello, World!";
        unsafe {
            std::ptr::copy_nonoverlapping(data.as_ptr(), mmap_file.fptr(), data.len());
        }
        mmap_file.load_file();
        mmap_file.close_file();
        remove_file("test.txt").unwrap();
    }

    #[test]
    fn test_sync_file() {
        let mut mmap_file = MMapFile::new();
        mmap_file.open_file("test.txt", 100).unwrap();
        let data = b"Hello, World!";
        unsafe {
            std::ptr::copy_nonoverlapping(data.as_ptr(), mmap_file.fptr(), data.len());
        }
        mmap_file.sync_file();
        mmap_file.close_file();

        let mut file = File::open("test.txt").unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        assert_eq!(contents, "Hello, World!");
        remove_file("test.txt").unwrap();
    }

    #[test]
    fn test_drop() {
        let file_path = "test.txt";
        {
            let mut mmap_file = MMapFile::new();
            mmap_file.open_file(file_path, 100).unwrap();
            let data = b"Hello, World!";
            unsafe {
                std::ptr::copy_nonoverlapping(data.as_ptr(), mmap_file.fptr(), data.len());
            }
        }
        assert!(File::open(file_path).is_ok());
        remove_file(file_path).unwrap();
    }
}
