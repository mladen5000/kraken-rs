use anyhow::{bail, Context, Result};
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}
impl Default for SequenceFormat {
    fn default() -> Self {
        SequenceFormat::AutoDetect
    }
}
#[derive(Default)]
pub struct Sequence {
    pub format: SequenceFormat,
    pub header: String,
    pub id: String,
    pub seq: String,
    pub quals: String,
    str_representation: String,
}

impl Sequence {
    pub fn to_string_repr(&mut self) -> &str {
        self.str_representation.clear();
        self.str_representation.push_str(&self.header);
        self.str_representation.push('\n');
        self.str_representation.push_str(&self.seq);
        self.str_representation.push('\n');
        if self.format == SequenceFormat::Fastq {
            self.str_representation.push_str("+\n");
            self.str_representation.push_str(&self.quals);
            self.str_representation.push('\n');
        }
        self.str_representation.as_ref()
    }
}

pub struct BatchSequenceReader<'a> {
    buffer: &'a [u8],
    buffer: Vec<String>, // Emulate reading chunks
}

impl BatchSequenceReader {
    pub fn new() -> Self {
        Self {
            file_format: SequenceFormat::AutoDetect,
            buffer: Vec::new(),
        }
    }

    fn strip_string(s: &mut String) {
        while s.ends_with(|c: char| c.is_whitespace()) {
            s.pop();
        }
    }

    // LoadBatch tries to read `record_count` FASTA or FASTQ records into internal buffer.
    // After detecting format, it reads either full records (for FASTQ: sets of 4 lines, for FASTA: until next '>').
    pub fn load_batch<R: BufRead>(&mut self, reader: &mut R, record_count: usize) -> Result<bool> {
        if self.file_format == SequenceFormat::AutoDetect {
            if let Some(&c) = reader.fill_buf()?.first() {
                let c = c as char;
                self.file_format = match c {
                    '@' => SequenceFormat::Fastq,
                    '>' => SequenceFormat::Fasta,
                    _ => bail!("sequence reader - unrecognized file format"),
                };
            } else {
                // empty input
                return Ok(false);
            }
        }

        self.buffer.clear();
        let mut records_left = record_count;
        let mut line_count = 0;
        let mut buf = String::new();
        while records_left > 0 {
            buf.clear();
            let n = reader.read_line(&mut buf)?;
            if n == 0 {
                break;
            }
            line_count += 1;
            self.buffer.push(buf.clone());
            match self.file_format {
                SequenceFormat::Fastq => {
                    // Every 4 lines is one record
                    if line_count % 4 == 0 {
                        records_left -= 1;
                    }
                }
                SequenceFormat::Fasta => {
                    // Every '>' indicates a new record
                    if let Some(&c) = reader.fill_buf()?.first() {
                        if c as char == '>' {
                            records_left -= 1;
                        }
                    } else {
                        // EOF
                        records_left = 0;
                    }
                }
                _ => {}
            }
        }

        Ok(!self.buffer.is_empty())
    }

    // Similar to LoadBlock logic is complicated due to streaming behavior.
    // We'll implement a simplified approach: read up to block_size bytes and then
    // try to detect format and store lines.
    pub fn load_block<R: Read>(&mut self, reader: &mut R, block_size: usize) -> Result<bool> {
        let mut block = vec![0u8; block_size];
        let count = reader.read(&mut block[..])?;
        if count == 0 {
            return Ok(false);
        }
        let content = String::from_utf8_lossy(&block[..count]).into_owned();
        let mut lines: Vec<_> = content.split('\n').map(|s| s.to_string()).collect();

        if self.file_format == SequenceFormat::AutoDetect {
            if let Some(first_line) = lines.get(0) {
                if let Some(c) = first_line.chars().next() {
                    self.file_format = match c {
                        '@' => SequenceFormat::Fastq,
                        '>' => SequenceFormat::Fasta,
                        _ => bail!("sequence reader - unrecognized file format"),
                    };
                } else {
                    return Ok(false);
                }
            } else {
                return Ok(false);
            }
        }

        // For FASTQ, we try to read through another few lines if possible (mimicking original)
        // For FASTA, we just store until next '>'.
        // The original code tries to top up lines after the block read.
        // We'll just store these lines for now.
        self.buffer = lines;
        Ok(!self.buffer.is_empty())
    }

    pub fn next_sequence(&mut self, seq: &mut Sequence) -> Result<bool> {
        Self::read_next_sequence(&mut self.buffer, seq, self.file_format)
    }

    fn read_next_sequence(
        lines: &mut Vec<String>,
        seq: &mut Sequence,
        mut file_format: SequenceFormat,
    ) -> Result<bool> {
        if lines.is_empty() {
            return Ok(false);
        }
        let mut line = lines.remove(0);
        Self::strip_string(&mut line);

        if file_format == SequenceFormat::AutoDetect {
            if let Some(c) = line.chars().next() {
                file_format = match c {
                    '@' => SequenceFormat::Fastq,
                    '>' => SequenceFormat::Fasta,
                    _ => bail!("sequence reader - unrecognized file format"),
                };
            } else {
                return Ok(false);
            }
        }

        seq.format = file_format;

        match seq.format {
            SequenceFormat::Fastq => {
                if line.is_empty() {
                    return Ok(false);
                }
                if !line.starts_with('@') {
                    bail!("malformed FASTQ file (expected '@', got '{}')", line);
                }
            }
            SequenceFormat::Fasta => {
                if !line.starts_with('>') {
                    bail!("malformed FASTA file (expected '>', got '{}')", line);
                }
            }
            _ => bail!("illegal sequence format encountered in parsing"),
        }

        seq.header = line.clone();
        let first_whitespace_ch = line[1..].find(|c: char| c.is_whitespace());
        // Adjust index because we did [1..]
        let substr_len = first_whitespace_ch.map(|idx| idx);

        if line.len() > 1 {
            if let Some(len) = substr_len {
                seq.id = line[1..(1 + len)].to_string();
            } else {
                // No whitespace found
                seq.id = line[1..].to_string();
            }
        } else {
            return Ok(false);
        }

        match seq.format {
            SequenceFormat::Fastq => {
                if lines.is_empty() {
                    return Ok(false);
                }
                let mut seq_line = lines.remove(0);
                Self::strip_string(&mut seq_line);
                seq.seq = seq_line;

                // Next line is '+', discard
                if lines.is_empty() {
                    return Ok(false);
                }
                let plus_line = lines.remove(0);
                if lines.is_empty() {
                    return Ok(false);
                }
                let mut quals_line = lines.remove(0);
                Self::strip_string(&mut quals_line);
                seq.quals = quals_line;
            }
            SequenceFormat::Fasta => {
                seq.quals.clear();
                seq.seq.clear();
                while !lines.is_empty() {
                    if let Some(peek_line) = lines.get(0) {
                        if peek_line.starts_with('>') {
                            // next sequence start
                            break;
                        }
                    }
                    let mut sline = lines.remove(0);
                    Self::strip_string(&mut sline);
                    seq.seq.push_str(&sline);
                }
            }
            _ => {}
        }

        Ok(true)
    }

    pub fn file_format(&self) -> SequenceFormat {
        self.file_format
    }
}

impl MinimizerScanner {
    // ...existing code...

    pub fn load_sequence(&mut self, seq: &str) {
        self.str_ = Some(seq.to_string());
        self.start_ = 0;
        self.finish_ = seq.len();
        self.str_pos_ = self.start_;
        if ((self.finish_ - self.start_) as isize) + 1 < self.l_ {
            // Invalidate scanner if interval < 1 l-mer
            self.str_pos_ = self.finish_;
        }
        self.queue_.clear();
        self.queue_pos_ = 0;
        self.loaded_ch_ = 0;
        self.last_minimizer_ = !0;
        self.last_ambig_ = 0;
        self.lmer_ = 0;
    }

    // ...existing code...
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_parsing() -> Result<()> {
        let data = "@read1\nACGT\n+\n!!!!\n@read2\nGATTACA\n+\n???????\n";
        let mut reader = BatchSequenceReader::new();
        {
            let mut bufreader = BufReader::new(data.as_bytes());
            reader.load_batch(&mut bufreader, 2)?;
        }
        let mut seq = Sequence {
            format: SequenceFormat::AutoDetect,
            header: String::new(),
            id: String::new(),
            seq: String::new(),
            quals: String::new(),
            str_representation: String::new(),
        };
        assert!(reader.next_sequence(&mut seq)?);
        assert_eq!(seq.id, "read1");
        assert_eq!(seq.seq, "ACGT");
        assert_eq!(seq.quals, "!!!!");
        assert!(reader.next_sequence(&mut seq)?);
        assert_eq!(seq.id, "read2");
        assert_eq!(seq.seq, "GATTACA");
        assert_eq!(seq.quals, "???????");
        assert!(!reader.next_sequence(&mut seq)?);
        Ok(())
    }

    #[test]
    fn test_fasta_parsing() -> Result<()> {
        let data = ">read1\nACGT\n>read2\nGATTACA\n";
        let mut reader = BatchSequenceReader::new();
        {
            let mut bufreader = BufReader::new(data.as_bytes());
            reader.load_batch(&mut bufreader, 2)?;
        }
        let mut seq = Sequence {
            format: SequenceFormat::AutoDetect,
            header: String::new(),
            id: String::new(),
            seq: String::new(),
            quals: String::new(),
            str_representation: String::new(),
        };
        assert!(reader.next_sequence(&mut seq)?);
        assert_eq!(seq.id, "read1");
        assert_eq!(seq.seq, "ACGT");
        assert!(reader.next_sequence(&mut seq)?);
        assert_eq!(seq.id, "read2");
        assert_eq!(seq.seq, "GATTACA");
        assert!(!reader.next_sequence(&mut seq)?);
        Ok(())
    }
}
