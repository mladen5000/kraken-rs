use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Removes trailing whitespace from a string.
pub fn strip_string(str: &mut String) {
    while str.ends_with(char::is_whitespace) {
        str.pop();
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}
/// Represents a biological sequence with its format, id, sequence data, and quality scores.
pub struct Sequence {
    pub format: SequenceFormat,
    pub id: String,
    pub seq: String,
    pub quals: String,
}

impl Sequence {
    pub fn new() -> Self {
        Sequence {
            format: SequenceFormat::AutoDetect,
            id: String::new(),
            seq: String::new(),
            quals: String::new(),
        }
    }

    /// Converts the sequence to a string representation based on its format.
    pub fn to_string(&self) -> String {
        let mut repr = self.id.clone();
        repr.push('\n');
        repr.push_str(&self.seq);
        match self.format {
            SequenceFormat::Fastq => {
                repr.push_str("\n+\n");
                repr.push_str(&self.quals);
            }
            _ => {}
        }
        repr.push('\n');
        repr
    }
}

/// Specifies the format of the sequence (FASTA, FASTQ, or auto-detection).

/// Reads sequences in batches from a file or stream.
pub struct BatchSequenceReader<R: BufRead> {
    reader: R,
    format: SequenceFormat,
}

impl<R: BufRead> BatchSequenceReader<R> {
    pub fn new(reader: R) -> Self {
        BatchSequenceReader {
            reader,
            format: SequenceFormat::AutoDetect,
        }
    }

    /// Loads a block of data from the input stream.
    pub fn load_block(&mut self, block_size: usize) -> io::Result<()> {
        let mut block = vec![0; block_size];
        let bytes_read = self.reader.read(&mut block)?;
        if bytes_read == 0 {
            return Ok(());
        }

        // Process block here...

        Ok(())
    }

    /// Loads a batch of sequences from the input stream.
    pub fn load_batch(&mut self, record_count: usize) -> io::Result<Vec<Sequence>> {
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
    pub fn read_next_sequence<B: BufRead>(
        reader: &mut R,
        file_format: &mut SequenceFormat,
    ) -> io::Result<Option<Sequence>> {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            return Ok(None); // End of input stream.
        }

        // Detect file format if set to AutoDetect.
        match file_format {
            SequenceFormat::AutoDetect => {
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
            _ => {}
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
            SequenceFormat::AutoDetect => todo!(),
        }

        Ok(Some(sequence))
    }

    /// Reads the next sequence from the stream.
    pub fn next_sequence(&mut self) -> io::Result<Option<Sequence>> {
        let mut line = String::new();

        // Detect file format if not already known.
        match self.format {
            SequenceFormat::AutoDetect => {
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
            }
            _ => {
                // Read the first line of the sequence (id).
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(None);
                }
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
