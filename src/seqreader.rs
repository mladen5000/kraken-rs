use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Removes trailing whitespace from a string.
fn strip_string(str: &mut String) {
    while str.ends_with(char::is_whitespace) {
        str.pop();
    }
}

/// Represents a biological sequence with its format, id, sequence data, and quality scores.
pub struct Sequence {
    pub format: SequenceFormat,
    pub id: String,
    pub seq: String,
    pub quals: String,
}

impl Sequence {
    fn new() -> Self {
        Sequence {
            format: SequenceFormat::AutoDetect,
            id: String::new(),
            seq: String::new(),
            quals: String::new(),
        }
    }

    /// Converts the sequence to a string representation based on its format.
    fn to_string(&self) -> String {
        let mut repr = self.id.clone();
        repr.push('\n');
        repr.push_str(&self.seq);
        if self.format == SequenceFormat::Fastq {
            repr.push_str("\n+\n");
            repr.push_str(&self.quals);
        }
        repr.push('\n');
        repr
    }
}

/// Specifies the format of the sequence (FASTA, FASTQ, or auto-detection).
pub enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}

/// Reads sequences in batches from a file or stream.
struct BatchSequenceReader<R: BufRead> {
    reader: R,
    format: SequenceFormat,
}

impl<R: BufRead> BatchSequenceReader<R> {
    fn new(reader: R) -> Self {
        BatchSequenceReader {
            reader,
            format: SequenceFormat::AutoDetect,
        }
    }

    /// Loads a block of data from the input stream.
    fn load_block(&mut self, block_size: usize) -> io::Result<()> {
        let mut block = vec![0; block_size];
        let bytes_read = self.reader.read(&mut block)?;
        if bytes_read == 0 {
            return Ok(());
        }

        // Process block here...

        Ok(())
    }

    /// Loads a batch of sequences from the input stream.
    fn load_batch(&mut self, record_count: usize) -> io::Result<Vec<Sequence>> {
        let mut sequences = Vec::new();

        for _ in 0..record_count {
            if let Some(sequence) = self.next_sequence()? {
                sequences.push(sequence);
            } else {
                break; // End of stream reached.
            }
        }

        Ok(sequences)
    }
    // Existing methods...

    /// Reads the next sequence from the given input stream.
    fn read_next_sequence<R: BufRead>(
        reader: &mut R,
        file_format: &mut SequenceFormat,
    ) -> io::Result<Option<Sequence>> {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            return Ok(None); // End of input stream.
        }

        // Detect file format if set to AutoDetect.
        if *file_format == SequenceFormat::AutoDetect {
            *file_format = match line.chars().next() {
                Some('@') => SequenceFormat::Fastq,
                Some('>') => SequenceFormat::Fasta,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Unrecognized file format",
                    ))
                }
            };
        }

        let mut sequence = Sequence::new();
        sequence.format = file_format.clone();
        sequence.id = line.trim_end().to_string();

        match file_format {
            SequenceFormat::Fasta => {
                sequence.seq.clear();
                while reader.read_line(&mut line)? != 0 {
                    if line.starts_with('>') {
                        break; // Start of next sequence.
                    }
                    sequence.seq.push_str(line.trim_end());
                }
            }
            SequenceFormat::Fastq => {
                // Read sequence line.
                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "Incomplete FASTQ record",
                    ));
                }
                sequence.seq = line.trim_end().to_string();

                // Skip plus line and read quality scores.
                line.clear();
                reader.read_line(&mut line)?;
                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "Incomplete FASTQ record",
                    ));
                }
                sequence.quals = line.trim_end().to_string();
            }
        }

        Ok(Some(sequence))
    }

    /// Reads the next sequence from the stream.
    fn next_sequence(&mut self) -> io::Result<Option<Sequence>> {
        let mut line = String::new();

        // Detect file format if not already known.
        if self.format == SequenceFormat::AutoDetect {
            if self.reader.read_line(&mut line)? == 0 {
                return Ok(None);
            }
            self.format = match line.chars().next() {
                Some('@') => SequenceFormat::Fastq,
                Some('>') => SequenceFormat::Fasta,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Unrecognized file format",
                    ))
                }
            };
        } else {
            // Read the first line of the sequence (id).
            if self.reader.read_line(&mut line)? == 0 {
                return Ok(None);
            }
        }
        strip_string(&mut line);

        let mut seq = Sequence::new();
        seq.format = self.format.clone();
        seq.id = line;

        match self.format {
            SequenceFormat::Fasta => {
                while self.reader.read_line(&mut line)? != 0 {
                    if line.starts_with('>') {
                        break; // Next sequence id reached.
                    }
                    strip_string(&mut line);
                    seq.seq.push_str(&line);
                    line.clear();
                }
            }
            SequenceFormat::Fastq => {
                // Read sequence line.
                if self.reader.read_line(&mut line)? == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "Incomplete FASTQ record",
                    ));
                }
                strip_string(&mut line);
                seq.seq = line;

                // Skip plus line.
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "Incomplete FASTQ record",
                    ));
                }

                // Read quality scores.
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "Incomplete FASTQ record",
                    ));
                }
                strip_string(&mut line);
                seq.quals = line;
            }
            SequenceFormat::AutoDetect => todo!(),
        }

        Ok(Some(seq))
    }
}

// Usage example:
fn main() -> io::Result<()> {
    let file = File::open("sequences.fasta")?;
    let reader = BufReader::new(file);
    let mut seq_reader = BatchSequenceReader::new(reader);

    while let Some(sequence) = seq_reader.next_sequence()? {
        println!("{}", sequence.to_string());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_fasta_sequence() {
        let data = ">seq1\nACGT\nTGCA\n";
        let mut reader = Cursor::new(data);
        let mut format = SequenceFormat::AutoDetect;

        let sequence = BatchSequenceReader::read_next_sequence(&mut reader, &mut format)
            .unwrap()
            .unwrap();
        assert_eq!(sequence.id, ">seq1");
        assert_eq!(sequence.seq, "ACGTTGCA");
        assert_eq!(sequence.quals, "");
        assert_eq!(format, SequenceFormat::Fasta);
    }

    #[test]
    fn test_read_fastq_sequence() {
        let data = "@seq1\nACGT\n+\n!!!!\n";
        let mut reader = Cursor::new(data);
        let mut format = SequenceFormat::AutoDetect;

        let sequence = BatchSequenceReader::read_next_sequence(&mut reader, &mut format)
            .unwrap()
            .unwrap();
        assert_eq!(sequence.id, "@seq1");
        assert_eq!(sequence.seq, "ACGT");
        assert_eq!(sequence.quals, "!!!!");
        assert_eq!(format, SequenceFormat::Fastq);
    }

    #[test]
    fn test_read_empty_sequence() {
        let data = "";
        let mut reader = Cursor::new(data);
        let mut format = SequenceFormat::AutoDetect;

        assert!(
            BatchSequenceReader::read_next_sequence(&mut reader, &mut format)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn test_invalid_format_detection() {
        let data = "invalid\nACGT\n";
        let mut reader = Cursor::new(data);
        let mut format = SequenceFormat::AutoDetect;

        assert!(BatchSequenceReader::read_next_sequence(&mut reader, &mut format).is_err());
    }

    // Additional tests for specific cases and error handling can be added here.
}
