use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;

pub struct GzMultiReader {
    files: Vec<String>,
    current_index: usize,
    current_reader: Option<GzDecoder<File>>,
    buffer: Vec<u8>,
    pos: usize,
    cap: usize,
}

impl GzMultiReader {
    pub fn new(filenames: Vec<String>) -> io::Result<Self> {
        let mut reader = GzMultiReader {
            files: filenames,
            current_index: 0,
            current_reader: None,
            buffer: vec![0; 8 * 1024], // 8KB buffer, similar to DEFAULT_BUFSIZ
            pos: 0,
            cap: 0,
        };
        reader.open_next_file()?;
        Ok(reader)
    }

    pub fn from_file<P: AsRef<Path>>(filename: P) -> io::Result<Self> {
        Self::new(vec![filename.as_ref().to_string_lossy().into_owned()])
    }

    fn open_next_file(&mut self) -> io::Result<()> {
        if self.current_index >= self.files.len() {
            self.current_reader = None;
            return Ok(());
        }

        let file = File::open(&self.files[self.current_index])?;
        self.current_reader = Some(GzDecoder::new(file));
        self.current_index += 1;
        self.pos = 0;
        self.cap = 0;
        Ok(())
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.cap {
            self.cap = match &mut self.current_reader {
                Some(reader) => reader.read(&mut self.buffer)?,
                None => return Ok(&[]),
            };
            self.pos = 0;
        }

        Ok(&self.buffer[self.pos..self.cap])
    }
}

impl Read for GzMultiReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut total_read = 0;

        while total_read < buf.len() {
            let available = self.fill_buf()?;
            if available.is_empty() {
                if self.open_next_file().is_ok() {
                    continue;
                } else {
                    break;
                }
            }

            let to_read = std::cmp::min(available.len(), buf.len() - total_read);
            buf[total_read..total_read + to_read].copy_from_slice(&available[..to_read]);
            self.pos += to_read;
            total_read += to_read;
        }

        Ok(total_read)
    }
}

impl BufRead for GzMultiReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.cap {
            self.cap = match &mut self.current_reader {
                Some(reader) => reader.read(&mut self.buffer)?,
                None => 0,
            };
            self.pos = 0;

            if self.cap == 0 && self.current_index < self.files.len() {
                self.open_next_file()?;
                return self.fill_buf();
            }
        }

        Ok(&self.buffer[self.pos..self.cap])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = std::cmp::min(self.pos + amt, self.cap);
    }
}
