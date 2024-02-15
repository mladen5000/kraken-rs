use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

pub struct GzInputStream {
    /// List of files to read from
    files: Vec<String>,
    /// Current file being read from
    current_file: Option<Box<dyn Read>>,
}

impl GzInputStream {
    /// Create a new GzInputStream from a single file
    pub fn new_single(filename: &str) -> GzInputStream {
        let file = File::open(filename).unwrap();
        let decoder = GzDecoder::new(file);
        let reader = Box::new(BufReader::new(decoder));

        GzInputStream {
            files: vec![filename.to_string()],
            current_file: Some(reader),
        }
    }

    /// Create a new GzInputStream from multiple files
    pub fn new_multiple(filenames: Vec<String>) -> GzInputStream {
        let file = File::open(&filenames[0]).unwrap();
        let decoder = GzDecoder::new(file);
        let reader = Box::new(BufReader::new(decoder));

        GzInputStream {
            files: filenames,
            current_file: Some(reader),
        }
    }

    pub fn next(&mut self) -> Option<u8> {
        match self.current_file.as_mut().unwrap().bytes().next() {
            Some(Ok(byte)) => Some(byte),
            _ => {
                self.files.remove(0);
                if !self.files.is_empty() {
                    let file = File::open(&self.files[0]).unwrap();
                    let decoder = GzDecoder::new(file);
                    self.current_file = Some(Box::new(BufReader::new(decoder)));
                    self.next()
                } else {
                    None
                }
            }
        }
    }
}
