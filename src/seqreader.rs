use std::fmt;
use std::io::{self, BufRead};

#[derive(Debug, Clone, Copy, PartialEq)]
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

fn strip_string(s: &mut String) {
    if s.is_empty() {
        return;
    }
    while s.ends_with(char::is_whitespace) {
        s.pop();
    }
}

#[derive(Default, Debug, Clone)]
pub struct Sequence {
    pub format: SequenceFormat,
    pub header: String, // header line, including @/>, but not newline
    pub id: String,     // from first char. after @/> up to first whitespace
    pub seq: String,
    pub quals: String, // only meaningful for FASTQ seqs
    pub str_representation: String,
}

impl Sequence {
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

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.format {
            SequenceFormat::Fastq => {
                writeln!(f, "{}", self.header)?;
                writeln!(f, "{}", self.seq)?;
                writeln!(f, "+")?;
                writeln!(f, "{}", self.quals)
            }
            _ => {
                writeln!(f, "{}", self.header)?;
                writeln!(f, "{}", self.seq)
            }
        }
    }
}

pub struct BatchSequenceReader<R: BufRead> {
    reader: R,
    str_buffer: String,
    block_buffer: Vec<u8>,
    file_format: SequenceFormat,
}

impl<R: BufRead> BatchSequenceReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            str_buffer: String::with_capacity(8192),
            block_buffer: vec![0; 8192],
            file_format: SequenceFormat::AutoDetect,
        }
    }

    pub fn file_format(&self) -> SequenceFormat {
        self.file_format
    }

    pub fn load_block(&mut self, block_size: usize) -> io::Result<bool> {
        if self.block_buffer.len() < block_size {
            self.block_buffer.resize(block_size, 0);
        }

        let bytes_read = self.reader.read(&mut self.block_buffer[..block_size])?;
        if bytes_read == 0 {
            return Ok(false);
        }

        // Auto-detect format if needed
        if self.file_format == SequenceFormat::AutoDetect {
            self.file_format = match self.block_buffer[0] as char {
                '@' => SequenceFormat::Fastq,
                '>' => SequenceFormat::Fasta,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "unrecognized file format",
                    ))
                }
            };
        }

        self.str_buffer.clear();
        self.str_buffer.push_str(
            std::str::from_utf8(&self.block_buffer[..bytes_read])
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        );

        // Read additional line to ensure we don't cut sequences
        let mut extra_line = String::new();
        if self.reader.read_line(&mut extra_line)? > 0 {
            self.str_buffer.push_str(&extra_line);
        }

        // For FASTQ, ensure we read complete records
        if self.file_format == SequenceFormat::Fastq {
            let mut lines_to_read = if extra_line.starts_with('@') { 3 } else { 2 };
            while lines_to_read > 0 && self.reader.read_line(&mut extra_line)? > 0 {
                self.str_buffer.push_str(&extra_line);
                extra_line.clear();
                lines_to_read -= 1;
            }
        } else {
            // For FASTA, read until next header or EOF
            loop {
                match self.reader.fill_buf()? {
                    buf if buf.is_empty() => break,
                    buf if buf[0] as char == '>' => break,
                    _ => {
                        extra_line.clear();
                        if self.reader.read_line(&mut extra_line)? == 0 {
                            break;
                        }
                        self.str_buffer.push_str(&extra_line);
                    }
                }
            }
        }

        Ok(true)
    }

    pub fn load_batch(&mut self, record_count: usize) -> io::Result<bool> {
        if record_count == 0 {
            return Ok(false);
        }

        self.str_buffer.clear();
        let mut valid = false;

        // Auto-detect format if needed
        if self.file_format == SequenceFormat::AutoDetect {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(false);
            }
            self.file_format = match buf[0] as char {
                '@' => SequenceFormat::Fastq,
                '>' => SequenceFormat::Fasta,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "unrecognized file format",
                    ))
                }
            };
            valid = true;
        }

        let mut line_count = 0;
        let mut records_read = 0;
        let mut line = String::new();

        while records_read < record_count {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                break;
            }
            line_count += 1;
            valid = true;

            self.str_buffer.push_str(&line);

            match self.file_format {
                SequenceFormat::Fastq if line_count % 4 == 0 => records_read += 1,
                SequenceFormat::Fasta => {
                    if self
                        .reader
                        .fill_buf()?
                        .first()
                        .map_or(false, |&c| c as char == '>')
                    {
                        records_read += 1;
                    }
                }
                _ => {}
            }
        }

        Ok(valid)
    }

    pub fn next_sequence(&mut self, sequence: &mut Sequence) -> io::Result<bool> {
        sequence.header.clear();
        sequence.id.clear();
        sequence.seq.clear();
        sequence.quals.clear();

        let mut line = String::new();
        if self.reader.read_line(&mut line)? == 0 {
            return Ok(false);
        }
        strip_string(&mut line);

        // Handle format detection and validation
        match (self.file_format, line.chars().next()) {
            (SequenceFormat::AutoDetect, Some('@')) => {
                self.file_format = SequenceFormat::Fastq;
                sequence.format = SequenceFormat::Fastq;
            }
            (SequenceFormat::AutoDetect, Some('>')) => {
                self.file_format = SequenceFormat::Fasta;
                sequence.format = SequenceFormat::Fasta;
            }
            (SequenceFormat::Fastq, Some('@')) => {
                sequence.format = SequenceFormat::Fastq;
            }
            (SequenceFormat::Fasta, Some('>')) => {
                sequence.format = SequenceFormat::Fasta;
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "malformed {} file",
                        if self.file_format == SequenceFormat::Fastq {
                            "FASTQ"
                        } else {
                            "FASTA"
                        }
                    ),
                ));
            }
        }

        sequence.header = line.clone();
        if sequence.header.len() > 1 {
            if let Some(end) = sequence.header[1..].find(|c: char| c.is_whitespace()) {
                sequence.id = sequence.header[1..end + 1].to_string();
            } else {
                sequence.id = sequence.header[1..].to_string();
            }
        } else {
            return Ok(false);
        }

        match sequence.format {
            SequenceFormat::Fastq => {
                // Read sequence
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(false);
                }
                strip_string(&mut line);
                sequence.seq = line.clone();

                // Read + line
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(false);
                }

                // Read quality scores
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(false);
                }
                strip_string(&mut line);
                sequence.quals = line.clone();
            }
            SequenceFormat::Fasta => {
                sequence.quals.clear();
                sequence.seq.clear();
                loop {
                    match self.reader.fill_buf()? {
                        buf if buf.is_empty() => break,
                        buf if buf[0] as char == '>' => break,
                        _ => {
                            line.clear();
                            if self.reader.read_line(&mut line)? == 0 {
                                break;
                            }
                            strip_string(&mut line);
                            sequence.seq.push_str(&line);
                        }
                    }
                }
            }
            _ => unreachable!(),
        }

        Ok(!sequence.seq.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parsing() {
        let data = ">read1 some desc\nACGT\n>read2 other desc\nGATTACA\n";
        let mut reader = BatchSequenceReader::new(data.as_bytes());
        let mut seq = Sequence::default();

        assert!(reader.next_sequence(&mut seq).unwrap());
        assert_eq!(seq.header, ">read1 some desc");
        assert_eq!(seq.id, "read1");
        assert_eq!(seq.seq, "ACGT");
        assert_eq!(seq.format, SequenceFormat::Fasta);

        assert!(reader.next_sequence(&mut seq).unwrap());
        assert_eq!(seq.header, ">read2 other desc");
        assert_eq!(seq.id, "read2");
        assert_eq!(seq.seq, "GATTACA");
        assert_eq!(seq.format, SequenceFormat::Fasta);

        assert!(!reader.next_sequence(&mut seq).unwrap());
    }

    #[test]
    fn test_fastq_parsing() {
        let data = "@read1 desc\nACGT\n+\n!!!!\n@read2 desc\nGATTACA\n+\n#######\n";
        let mut reader = BatchSequenceReader::new(data.as_bytes());
        let mut seq = Sequence::default();

        assert!(reader.next_sequence(&mut seq).unwrap());
        assert_eq!(seq.header, "@read1 desc");
        assert_eq!(seq.id, "read1");
        assert_eq!(seq.seq, "ACGT");
        assert_eq!(seq.quals, "!!!!");
        assert_eq!(seq.format, SequenceFormat::Fastq);

        assert!(reader.next_sequence(&mut seq).unwrap());
        assert_eq!(seq.header, "@read2 desc");
        assert_eq!(seq.id, "read2");
        assert_eq!(seq.seq, "GATTACA");
        assert_eq!(seq.quals, "#######");
        assert_eq!(seq.format, SequenceFormat::Fastq);

        assert!(!reader.next_sequence(&mut seq).unwrap());
    }
}
