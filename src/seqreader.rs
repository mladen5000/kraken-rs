/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::io::{BufRead, Cursor, Read};
use std::str::from_utf8;

#[derive(PartialEq, PartialOrd, Copy, Clone, Debug)]
pub enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct Sequence {
    pub format: SequenceFormat,
    pub header: String,
    pub id: String,
    pub seq: String,
    pub quals: String,
    str_representation: String,
}

impl Default for Sequence {
    fn default() -> Self {
        Self::new()
    }
}

impl Sequence {
    pub fn new() -> Self {
        Sequence {
            format: SequenceFormat::AutoDetect,
            header: "".to_string(),
            id: "".to_string(),
            seq: "".to_string(),
            quals: "".to_string(),
            str_representation: "".to_string(),
        }
    }
    pub fn to_string(&mut self) -> &str {
        self.str_representation.clear();
        self.str_representation.push_str(&self.header);

        match self.format {
            SequenceFormat::Fastq => {
                self.str_representation.push('\n');
                self.str_representation.push_str(&self.seq);
                self.str_representation.push_str("\n+\n");
                self.str_representation.push_str(&self.quals);
                self.str_representation.push('\n');
            }
            _ => {
                self.str_representation.push('\n');
                self.str_representation.push_str(&self.seq);
                self.str_representation.push('\n');
            }
        }

        &self.str_representation
    }
}

pub struct BatchSequenceReader {
    // Sequence buffer
    ss: Vec<u8>,
    // String buffer
    str_buffer: String,
    // File format
    file_format: SequenceFormat,
    // Block buffer
    block_buffer: Vec<u8>,
}

impl Default for BatchSequenceReader {
    fn default() -> Self {
        Self::new()
    }
}

impl BatchSequenceReader {
    pub fn new() -> Self {
        BatchSequenceReader {
            ss: Vec::with_capacity(8192),
            str_buffer: String::with_capacity(8192),
            file_format: SequenceFormat::AutoDetect,
            block_buffer: vec![0; 8192],
        }
    }

    pub fn load_block<R: BufRead>(&mut self, ifs: &mut R, block_size: usize) -> bool {
        self.ss.clear();

        if self.block_buffer.len() < block_size {
            self.block_buffer.resize(block_size, 0);
        }

        let bytes_read = ifs.read(&mut self.block_buffer[..block_size]).unwrap_or(0);
        if bytes_read == 0 {
            return false;
        }

        if self.file_format == SequenceFormat::AutoDetect {
            match self.block_buffer[0] {
                b'@' => self.file_format = SequenceFormat::Fastq,
                b'>' => self.file_format = SequenceFormat::Fasta,
                _ => panic!("sequence reader - unrecognized file format"),
            }
        }

        self.str_buffer.clear();
        self.str_buffer
            .push_str(from_utf8(&self.block_buffer[..bytes_read]).unwrap());
        self.ss.extend_from_slice(self.str_buffer.as_bytes());

        if let Some(line) = ifs.lines().next() {
            self.str_buffer.clear();
            self.str_buffer.push_str(&line.unwrap());
            self.str_buffer.push('\n');
            self.ss.extend_from_slice(self.str_buffer.as_bytes());
        }

        if self.file_format == SequenceFormat::Fastq {
            while let Some(line) = ifs.lines().next() {
                self.str_buffer.clear();
                self.str_buffer.push_str(&line.unwrap());
                self.str_buffer.push('\n');
                self.ss.extend_from_slice(self.str_buffer.as_bytes());
                if self.str_buffer.as_bytes()[0] == b'@' {
                    break;
                }
            }

            let lines_to_read;
            if let Some(line) = ifs.lines().next() {
                self.str_buffer.clear();
                self.str_buffer.push_str(&line.unwrap());
                self.str_buffer.push('\n');
                self.ss.extend_from_slice(self.str_buffer.as_bytes());
                lines_to_read = if self.str_buffer.as_bytes()[0] == b'@' {
                    3
                } else {
                    2
                };
                for _ in 0..lines_to_read {
                    if let Some(line) = ifs.lines().next() {
                        self.str_buffer.clear();
                        self.str_buffer.push_str(&line.unwrap());
                        self.str_buffer.push('\n');
                        self.ss.extend_from_slice(self.str_buffer.as_bytes());
                    }
                }
            }
        } else {
            while let Some(byte) = ifs.bytes().next() {
                if byte.unwrap() == b'>' {
                    break;
                }
                if let Some(line) = ifs.lines().next() {
                    self.str_buffer.clear();
                    self.str_buffer.push_str(&line.unwrap());
                    self.str_buffer.push('\n');
                    self.ss.extend_from_slice(self.str_buffer.as_bytes());
                }
            }
        }

        true
    }

    pub fn load_batch<R: BufRead>(&mut self, ifs: &mut R, record_count: usize) -> bool {
        self.ss.clear();

        let mut valid = false;
        if self.file_format == SequenceFormat::AutoDetect {
            if let Some(byte) = ifs.bytes().next() {
                match byte.unwrap() {
                    b'@' => self.file_format = SequenceFormat::Fastq,
                    b'>' => self.file_format = SequenceFormat::Fasta,
                    _ => panic!("sequence reader - unrecognized file format"),
                }
                valid = true;
            } else {
                return false;
            }
        }

        let mut line_count = 0;
        let mut remaining_record_count = record_count;
        while remaining_record_count > 0 {
            if let Some(line) = ifs.lines().next() {
                self.str_buffer.clear();
                self.str_buffer.push_str(&line.unwrap());
                self.str_buffer.push('\n');
                self.ss.extend_from_slice(self.str_buffer.as_bytes());
                line_count += 1;
                valid = true;

                if self.file_format == SequenceFormat::Fastq {
                    if line_count % 4 == 0 {
                        remaining_record_count -= 1;
                    }
                } else if let Some(byte) = ifs.bytes().next() {
                    if byte.unwrap() == b'>' {
                        remaining_record_count -= 1;
                    }
                }
            } else {
                break;
            }
        }

        valid
    }

    pub fn next_sequence(&mut self, seq: &mut Sequence) -> Option<bool> {
        let mut cursor = Cursor::new(&mut self.ss);
        let result = BatchSequenceReader::read_next_sequence(
            &mut cursor, // input stream
            seq,
            &mut self.str_buffer,
            self.file_format,
        );

        if result {
            let pos = cursor.position() as usize;
            self.ss.drain(..pos);
            Some(result)
        } else {
            None
        }
    }

    /// Read the next sequence from the input stream
    pub fn read_next_sequence<R: BufRead>(
        is: &mut R, // input stream
        seq: &mut Sequence,
        str_buffer: &mut String,
        file_format: SequenceFormat,
    ) -> bool {
        // str_buffer.clear();
        println!("{}", str_buffer.len());

        if is.read_line(str_buffer).unwrap() == 0 {
            return false;
        }
        strip_string(str_buffer);

        let mut format = file_format;
        if format == SequenceFormat::AutoDetect {
            match str_buffer.as_bytes()[0] {
                b'@' => format = SequenceFormat::Fastq,
                b'>' => format = SequenceFormat::Fasta,
                _ => panic!("sequence reader - unrecognized file format"),
            }
        }

        seq.format = format;
        if seq.format == SequenceFormat::Fastq {
            if str_buffer.is_empty() {
                return false;
            }
            if str_buffer.as_bytes()[0] != b'@' {
                panic!(
                    "malformed FASTQ file (exp. '@', saw \"{}\"), aborting",
                    str_buffer
                );
            }
        } else if seq.format == SequenceFormat::Fasta {
            if str_buffer.as_bytes()[0] != b'>' {
                panic!(
                    "malformed FASTA file (exp. '>', saw \"{}\"), aborting",
                    str_buffer
                );
            }
        } else {
            panic!("illegal sequence format encountered in parsing");
        }

        seq.header.clear();
        seq.header.push_str(str_buffer);

        let first_whitespace_ch = str_buffer
            .find(|c: char| c.is_whitespace())
            .unwrap_or(str_buffer.len());
        let substr_len = if first_whitespace_ch > 1 {
            first_whitespace_ch // - 1
        } else {
            str_buffer.len() - 1
        };

        if substr_len > 0 {
            seq.id.clear();
            seq.id.push_str(&str_buffer[1..=substr_len]);
        } else {
            seq.id.clear();
        }

        seq.seq.clear();
        if seq.format == SequenceFormat::Fasta {
            while let Ok(bytes_read) = is.read_line(str_buffer) {
                if bytes_read == 0 || str_buffer.starts_with('>') {
                    break;
                }
                strip_string(str_buffer);
                seq.seq.push_str(str_buffer);
                str_buffer.clear();
            }
        } else if seq.format == SequenceFormat::Fastq {
            if is.read_line(str_buffer).unwrap() == 0 {
                return false;
            }
            strip_string(str_buffer);
            seq.seq.push_str(str_buffer);

            if is.read_line(str_buffer).unwrap() == 0 {
                return false;
            }

            if is.read_line(str_buffer).unwrap() == 0 {
                return false;
            }
            strip_string(str_buffer);
            seq.quals.clear();
            seq.quals.push_str(str_buffer);
        }

        true
    }
    pub fn file_format(&self) -> SequenceFormat {
        self.file_format
    }
}

fn strip_string(str: &mut String) {
    if str.is_empty() {
        return;
    }
    while let Some(c) = str.pop() {
        if !c.is_whitespace() {
            str.push(c);
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_sequence_to_string() {
        let mut seq = Sequence {
            format: SequenceFormat::Fasta,
            header: ">seq1".to_string(),
            id: "seq1".to_string(),
            seq: "ATCG".to_string(),
            quals: "".to_string(),
            str_representation: "".to_string(),
        };
        assert_eq!(seq.to_string(), ">seq1\nATCG\n");

        seq.format = SequenceFormat::Fastq;
        seq.quals = "!@#$".to_string();
        assert_eq!(seq.to_string(), ">seq1\nATCG\n+\n!@#$\n");
    }

    #[test]
    fn test_load_block_fasta() {
        let fasta_data = b">seq1\nATCG\n>seq2\nTGCA\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fasta_data);
        assert!(reader.load_block(&mut cursor, fasta_data.len()));
        assert_eq!(reader.file_format(), SequenceFormat::Fasta);
    }

    #[test]
    fn test_load_block_fastq() {
        let fastq_data = b"@seq1\nATCG\n+\n!@#$\n@seq2\nTGCA\n+\n$#@!\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fastq_data);
        assert!(reader.load_block(&mut cursor, fastq_data.len()));
        assert_eq!(reader.file_format(), SequenceFormat::Fastq);
    }

    #[test]
    fn test_load_batch_fasta() {
        let fasta_data = b">seq1\nATCG\n>seq2\nTGCA\n>seq3\nGCAT\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fasta_data);
        assert!(reader.load_batch(&mut cursor, 2));
        assert_eq!(reader.file_format(), SequenceFormat::Fasta);
    }

    #[test]
    fn test_load_batch_fastq() {
        let fastq_data = b"@seq1\nATCG\n+\n!@#$\n@seq2\nTGCA\n+\n$#@!\n@seq3\nGCAT\n+\n!@#$\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fastq_data);
        assert!(reader.load_batch(&mut cursor, 2));
        assert_eq!(reader.file_format(), SequenceFormat::Fastq);
    }

    #[test]
    fn test_next_sequence_fasta() {
        let fasta_data = b">seq1\nATCG\n>seq2\nTGCA\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fasta_data);
        reader.load_block(&mut cursor, fasta_data.len());

        let mut seq = Sequence {
            format: SequenceFormat::AutoDetect,
            header: "".to_string(),
            id: "".to_string(),
            seq: "".to_string(),
            quals: "".to_string(),
            str_representation: "".to_string(),
        };

        assert!(reader.next_sequence(&mut seq).is_some());
        assert_eq!(seq.format, SequenceFormat::Fasta);
        assert_eq!(seq.header, ">seq1");
        assert_eq!(seq.id, "seq1");
        assert_eq!(seq.seq, "ATCG");
        assert_eq!(seq.quals, "");

        assert!(reader.next_sequence(&mut seq).is_some());
        assert_eq!(seq.format, SequenceFormat::Fasta);
        assert_eq!(seq.header, ">seq2");
        assert_eq!(seq.id, "seq2");
        assert_eq!(seq.seq, "TGCA");
        assert_eq!(seq.quals, "");

        assert!(reader.next_sequence(&mut seq).is_none());
    }

    #[test]
    fn test_next_sequence_fastq() {
        let fastq_data = b"@seq1\nATCG\n+\n!@#$\n@seq2\nTGCA\n+\n$#@!\n";
        let mut reader = BatchSequenceReader::new();
        let mut cursor = Cursor::new(fastq_data);
        reader.load_block(&mut cursor, fastq_data.len());

        let mut seq = Sequence {
            format: SequenceFormat::AutoDetect,
            header: "".to_string(),
            id: "".to_string(),
            seq: "".to_string(),
            quals: "".to_string(),
            str_representation: "".to_string(),
        };

        assert!(reader.next_sequence(&mut seq).is_some());
        assert_eq!(seq.format, SequenceFormat::Fastq);
        assert_eq!(seq.header, "@seq1");
        assert_eq!(seq.id, "seq1");
        assert_eq!(seq.seq, "ATCG");
        assert_eq!(seq.quals, "!@#$");

        assert!(reader.next_sequence(&mut seq).is_some());
        assert_eq!(seq.format, SequenceFormat::Fastq);
        assert_eq!(seq.header, "@seq2");
        assert_eq!(seq.id, "seq2");
        assert_eq!(seq.seq, "TGCA");
        assert_eq!(seq.quals, "$#@!");
    }

    #[test]
    fn test_strip_string() {
        let mut s1 = "  hello  ".to_string();
        strip_string(&mut s1);
        assert_eq!(s1, "  hello");

        let mut s2 = "world  ".to_string();
        strip_string(&mut s2);
        assert_eq!(s2, "world");

        let mut s3 = "".to_string();
        strip_string(&mut s3);
        assert_eq!(s3, "");
    }
}
