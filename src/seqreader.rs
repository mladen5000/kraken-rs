use std::io::{BufRead, BufReader, Read};

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum SequenceFormat {
    AutoDetect,
    Fasta,
    Fastq,
}

#[derive(Clone)]
pub struct Sequence {
    pub format: SequenceFormat,
    pub header: String, // Header line, including @/>, but not newline
    pub id: String,     // From first char after @/> up to first whitespace
    pub seq: String,
    pub quals: String, // Only meaningful for FASTQ sequences
    pub str_representation: String,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            format: SequenceFormat::AutoDetect,
            header: String::new(),
            id: String::new(),
            seq: String::new(),
            quals: String::new(),
            str_representation: String::new(),
        }
    }

    pub fn to_string(&mut self) -> &String {
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

pub struct BatchSequenceReader<R: Read> {
    reader: BufReader<R>,
    block_size: usize,
    buffer: String,
    file_format: SequenceFormat,
    block_buffer: Vec<u8>,
}

impl<R: Read> BatchSequenceReader<R> {
    pub fn new(reader: R, block_size: usize) -> Self {
        Self {
            reader: BufReader::new(reader),
            block_size,
            buffer: String::with_capacity(block_size),
            file_format: SequenceFormat::AutoDetect,
            block_buffer: vec![0; block_size],
        }
    }

    pub fn load_block(&mut self, reader: &mut dyn Read, block_size: usize) -> bool {
        self.buffer.clear();
        if block_size > self.block_buffer.len() {
            self.block_buffer.resize(block_size, 0);
        }

        match reader.read(&mut self.block_buffer[..block_size]) {
            Ok(0) => return false,
            Ok(bytes_read) => {
                self.buffer
                    .push_str(&String::from_utf8_lossy(&self.block_buffer[..bytes_read]));
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => return false,
            Err(e) => panic!("Failed to read block: {:?}", e),
        }

        if self.file_format == SequenceFormat::AutoDetect {
            self.file_format = match self.block_buffer[0] {
                b'@' => SequenceFormat::Fastq,
                b'>' => SequenceFormat::Fasta,
                _ => panic!("Sequence reader - unrecognized file format"),
            };
        }

        // Convert &str to &[u8] and pass it to BufReader
        let buf = self.buffer.clone(); // Clone buffer to avoid mutable borrowing
        let mut buf_reader = BufReader::new(buf.as_bytes()); // Now use the clone's as_bytes
        if self.file_format == SequenceFormat::Fastq {
            self.read_fastq_sequences(&mut buf_reader)
        } else {
            self.read_fasta_sequences(&mut buf_reader)
        }
    }

    pub fn next_sequence(&mut self, seq: &mut Sequence) -> bool {
        let buffer_clone = self.buffer.clone(); // Clone the buffer to avoid mutable borrowing issues
        let format = self.file_format.clone(); // Clone the format to avoid moving it

        // Create a BufReader from the cloned buffer slice
        let mut buf_reader = BufReader::new(buffer_clone.as_bytes());

        BatchSequenceReader::read_next_sequence(
            &mut buf_reader, // Use the clone for the buffer
            seq,
            format,
        )
    }

    fn read_next_sequence(
        reader: &mut BufReader<&[u8]>,
        sequence: &mut Sequence,
        format: SequenceFormat,
    ) -> bool {
        let mut header = String::new();

        // Read the header line
        if reader.read_line(&mut header).unwrap() == 0 {
            return false; // EOF or empty line
        }

        // Remove trailing newline or carriage return characters
        header = header.trim_end().to_string();

        // Split header into ID and remaining header info
        if format == SequenceFormat::Fasta {
            if header.starts_with('>') {
                println!("header: {:?}", header);
                sequence.id = header[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
            } else {
                eprintln!("Error: Expected '>' at the beginning of the FASTA header.");
                return false;
            }
        } else if format == SequenceFormat::Fastq {
            if header.starts_with('@') {
                sequence.id = header[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
            } else {
                eprintln!("Error: Expected '@' at the beginning of the FASTQ header.");
                return false;
            }
        }

        // Read the sequence lines
        let mut sequence_data = String::new();
        let mut line = String::new();

        while reader.read_line(&mut line).unwrap() > 0 {
            line = line.trim_end().to_string();

            if line.starts_with('>') || line.starts_with('@') || line.starts_with('+') {
                // This line indicates the start of the next sequence or a section in FASTQ
                // Save the sequence data and return true, indicating that a sequence was successfully read.
                sequence.seq = sequence_data;
                return true; // Indicate that a sequence has been successfully read
            }

            sequence_data.push_str(&line);
            line = String::new(); // Reset the line for the next iteration
        }

        // If we reach here, it means we've finished reading the file without encountering a new sequence
        if !sequence_data.is_empty() {
            sequence.seq = sequence_data;
            return true; // Return true as the last sequence has been read
        }

        sequence.seq = sequence_data;

        // If FASTQ, skip the '+' line and read the quality line
        if format == SequenceFormat::Fastq {
            let mut plus_line = String::new();
            reader.read_line(&mut plus_line).unwrap(); // Read '+' line

            let mut quals = String::new();
            reader.read_line(&mut quals).unwrap(); // Read quality line
            sequence.quals = quals.trim_end().to_string();
        }

        true
    }

    /*
    fn some_function() {
        let buffer = b">seq1\nATCGATCGATCG\n>seq2\nGCTAGCTAGCTA\n";
        let mut buf_reader = BufReader::new(&buffer[..]); // Create a BufReader from the byte slice
        let mut sequence = Sequence::default();
        let format = SequenceFormat::Fasta;

        let _ = BatchSequenceReader::read_next_sequence(&mut buf_reader, &mut sequence, format);
    }
    */

    fn read_fastq_sequences(&mut self, reader: &mut BufReader<&[u8]>) -> bool {
        while let Ok(_) = reader.read_line(&mut self.buffer) {
            if self.buffer.is_empty() {
                break;
            }
            self.buffer.push_str("\n");
            if self.buffer.starts_with('@') {
                break;
            }
        }
        true
    }

    fn read_fasta_sequences(&mut self, reader: &mut BufReader<&[u8]>) -> bool {
        while let Ok(_) = reader.read_line(&mut self.buffer) {
            if self.buffer.is_empty() {
                break;
            }
            self.buffer.push_str("\n");
            if self.buffer.starts_with('>') {
                break;
            }
        }
        true
    }

    pub fn file_format(&self) -> SequenceFormat {
        self.file_format.clone() // Return a clone of the format to avoid moving it
    }
}

fn strip_string(s: &mut String) {
    while s.ends_with(char::is_whitespace) {
        s.pop();
    }
}
pub fn main() {
    // Example sequence in FASTA format
    let test_sequence = ">TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n";

    // Initialize the sequence and reader
    let mut sequence = Sequence::new();
    let mut reader = BatchSequenceReader::new();

    // Simulate reading the sequence from a buffer
    reader.buffer.push_str(test_sequence);

    // Load the sequence and process it
    if reader.next_sequence(&mut sequence) {
        println!("Processed Sequence:");
        println!("{}", sequence.to_string());
    } else {
        println!("Failed to process sequence.");
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_sequence() {
        let fasta_data = ">TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n";
        let mut sequence = Sequence::new();
        let mut reader = BatchSequenceReader::new();

        reader.buffer.push_str(fasta_data);

        assert!(
            reader.next_sequence(&mut sequence),
            "Failed to read sequence."
        );
        assert_eq!(
            sequence.id, "TestSequence",
            "ID mismatch: {:?}",
            sequence.id
        );
        assert_eq!(
            sequence.seq, "AGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT",
            "Sequence mismatch: {:?}",
            sequence.seq
        );
        assert_eq!(
            sequence.format,
            SequenceFormat::Fasta,
            "Format mismatch: {:?}",
            sequence.format
        );
    }

    #[test]
    fn test_fastq_sequence() {
        let fastq_data = "@TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        let mut sequence = Sequence::new();
        let mut reader = BatchSequenceReader::new();

        reader.buffer.push_str(fastq_data);

        assert!(
            reader.next_sequence(&mut sequence),
            "Failed to read sequence."
        );
        assert_eq!(
            sequence.id, "TestSequence",
            "ID mismatch: {:?}",
            sequence.id
        );
        assert_eq!(
            sequence.seq, "AGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT",
            "Sequence mismatch: {:?}",
            sequence.seq
        );
        assert_eq!(
            sequence.quals, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "Quals mismatch: {:?}",
            sequence.quals
        );
        assert_eq!(
            sequence.format,
            SequenceFormat::Fastq,
            "Format mismatch: {:?}",
            sequence.format
        );
    }

    #[test]
    fn test_fasta_to_string() {
        let fasta_data = ">TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n";
        let mut sequence = Sequence::new();
        let mut reader = BatchSequenceReader::new();

        reader.buffer.push_str(fasta_data);

        assert!(
            reader.next_sequence(&mut sequence),
            "Failed to read sequence."
        );
        let expected_output = ">TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n";
        assert_eq!(
            sequence.to_string().clone(),
            expected_output,
            "Output mismatch: {:?}",
            sequence.to_string()
        );
    }

    #[test]
    fn test_fastq_to_string() {
        let fastq_data = "@TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        let mut sequence = Sequence::new();
        let mut reader = BatchSequenceReader::new();

        reader.buffer.push_str(fastq_data);

        assert!(
            reader.next_sequence(&mut sequence),
            "Failed to read sequence."
        );
        let expected_output = "@TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        assert_eq!(
            sequence.to_string().clone(),
            expected_output,
            "Output mismatch: {:?}",
            sequence.to_string()
        );
    }

    #[test]
    fn test_malformed_fastq() {
        let malformed_fastq_data = "@TestSequence\nAGCTGATCGTAGCTAGCTAGCTGATCGTAGCTGACT\n+\n";
        let mut sequence = Sequence::new();
        let mut reader = BatchSequenceReader::new();

        reader.buffer.push_str(malformed_fastq_data);

        assert!(
            reader.next_sequence(&mut sequence),
            "Failed to read sequence."
        );
        assert!(
            sequence.quals.is_empty(),
            "Quals should be empty for malformed FASTQ, got: {:?}",
            sequence.quals
        );
    }
}
