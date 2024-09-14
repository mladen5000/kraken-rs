/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::fs::OpenOptions;
use std::io::{self, ErrorKind};
use std::mem;
use std::os::unix::fs::OpenOptionsExt;
use std::os::unix::io::AsRawFd;
use std::path::Path;
use std::ptr;
use std::slice;

use libc::{
    c_void, close, fstat, lseek, mmap, munmap, off_t, size_t, stat, sysconf, write as libc_write,
    MAP_FAILED, MAP_PRIVATE, MAP_SHARED, MS_SYNC, O_ACCMODE, O_APPEND, O_CREAT, O_RDONLY, O_RDWR,
    O_WRONLY, PROT_READ, PROT_WRITE, SEEK_SET, _SC_PAGESIZE,
};

pub struct MMapFile {
    ptr: *mut u8,
    size: usize,
    fd: i32,
    valid: bool,
}

impl MMapFile {
    /// Creates a new, uninitialized `MMapFile`.
    pub fn new() -> Self {
        MMapFile {
            ptr: ptr::null_mut(),
            size: 0,
            fd: -1,
            valid: false,
        }
    }

    /// Opens a file and memory-maps it.
    ///
    /// # Arguments
    ///
    /// * `filename` - The path to the file.
    /// * `mode` - The file open mode flags (e.g., `O_RDONLY`, `O_RDWR`).
    /// * `map_flags` - Optional mapping flags (e.g., `MAP_PRIVATE`, `MAP_SHARED`).
    /// * `prot_flags` - Optional protection flags (e.g., `PROT_READ`, `PROT_WRITE`).
    /// * `size` - Optional size of the mapping. Required if `O_CREAT` is set in `mode`.
    ///
    /// # Returns
    ///
    /// * `io::Result<()>` indicating success or failure.
    pub fn open<P: AsRef<Path>>(
        &mut self,
        filename: P,
        mode: i32,
        map_flags: Option<i32>,
        prot_flags: Option<i32>,
        size: Option<usize>,
    ) -> io::Result<()> {
        let filename = filename.as_ref();
        if (mode & O_APPEND != 0) || (mode & O_ACCMODE == O_WRONLY) {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "Illegal mode passed to MMapFile",
            ));
        }

        let prot_flags = prot_flags.unwrap_or_else(|| {
            if mode & O_ACCMODE == O_RDONLY {
                PROT_READ
            } else {
                PROT_READ | PROT_WRITE
            }
        });

        let map_flags = map_flags.unwrap_or_else(|| {
            if mode & O_ACCMODE == O_RDONLY {
                MAP_PRIVATE
            } else {
                MAP_SHARED
            }
        });

        let file = OpenOptions::new()
            .read(mode & O_ACCMODE != O_WRONLY)
            .write(mode & O_ACCMODE != O_RDONLY)
            .custom_flags(mode & !O_ACCMODE)
            .open(filename)?;

        self.fd = file.as_raw_fd();

        if mode & O_CREAT != 0 {
            let size = size.ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    "Size must be provided when creating a new file",
                )
            })?;
            unsafe {
                if lseek(self.fd, (size - 1) as off_t, SEEK_SET) < 0 {
                    return Err(io::Error::last_os_error());
                }
                // Write a single null byte to set the file size.
                if libc_write(self.fd, b"\0".as_ptr() as *const c_void, 1) < 0 {
                    return Err(io::Error::last_os_error());
                }
            }
            self.size = size;
        } else {
            let mut stat_buf: stat = unsafe { mem::zeroed() };
            let res = unsafe { fstat(self.fd, &mut stat_buf) };
            if res < 0 {
                return Err(io::Error::last_os_error());
            }
            self.size = stat_buf.st_size as usize;
        }

        self.ptr = unsafe {
            mmap(
                ptr::null_mut(),
                self.size as size_t,
                prot_flags,
                map_flags,
                self.fd,
                0,
            ) as *mut u8
        };

        if self.ptr == MAP_FAILED as *mut u8 {
            return Err(io::Error::last_os_error());
        }

        self.valid = true;
        Ok(())
    }

    /// Provides immutable access to the memory-mapped data as a byte slice.
    ///
    /// # Returns
    ///
    /// * `Option<&[u8]>` containing the data slice if valid, or `None` otherwise.
    pub fn as_slice(&self) -> Option<&[u8]> {
        if self.valid {
            unsafe { Some(slice::from_raw_parts(self.ptr, self.size)) }
        } else {
            None
        }
    }

    /// Returns the size of the memory-mapped file.
    ///
    /// # Returns
    ///
    /// * `usize` representing the file size in bytes.
    pub fn filesize(&self) -> usize {
        if self.valid {
            self.size
        } else {
            0
        }
    }

    /// Loads the file into the OS cache by reading it in pages.
    ///
    /// # Returns
    ///
    /// * `io::Result<()>` indicating success or failure.
    pub fn load_file(&self) -> io::Result<()> {
        if !self.valid {
            return Err(io::Error::new(
                ErrorKind::Other,
                "Cannot load an invalid mmap file",
            ));
        }

        let page_size = unsafe { sysconf(_SC_PAGESIZE) as usize };
        let num_pages = (self.size + page_size - 1) / page_size;

        let mut buf = vec![0u8; page_size];

        for i in 0..num_pages {
            let offset = i * page_size;
            let size = if offset + page_size <= self.size {
                page_size
            } else {
                self.size - offset
            };
            unsafe {
                ptr::copy_nonoverlapping(self.ptr.add(offset), buf.as_mut_ptr(), size);
            }
        }

        Ok(())
    }

    /// Synchronizes the memory-mapped file with the underlying storage.
    ///
    /// # Returns
    ///
    /// * `io::Result<()>` indicating success or failure.
    pub fn sync_file(&self) -> io::Result<()> {
        if self.valid {
            let res = unsafe { libc::msync(self.ptr as *mut c_void, self.size, MS_SYNC) };
            if res != 0 {
                return Err(io::Error::last_os_error());
            }
            Ok(())
        } else {
            Err(io::Error::new(
                ErrorKind::Other,
                "Cannot sync an invalid mmap file",
            ))
        }
    }

    /// Closes the memory-mapped file, ensuring data is synchronized and resources are freed.
    ///
    /// # Returns
    ///
    /// * `io::Result<()>` indicating success or failure.
    pub fn close_file(&mut self) -> io::Result<()> {
        if self.valid {
            self.sync_file()?;
            unsafe {
                munmap(self.ptr as *mut c_void, self.size);
                close(self.fd);
            }
            self.valid = false;
            Ok(())
        } else {
            Ok(())
        }
    }
}

impl Drop for MMapFile {
    fn drop(&mut self) {
        let _ = self.close_file();
    }
}
