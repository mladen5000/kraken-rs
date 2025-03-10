// Constants and helper function.
pub const KS_SEP_SPACE: i32 = 0; // whitespace (e.g. '\t', '\n', etc.)
pub const KS_SEP_TAB: i32 = 1; // whitespace except for ' '
pub const KS_SEP_LINE: i32 = 2; // newline (Unix/Windows)
pub const KS_SEP_MAX: i32 = 2;

pub fn kroundup32(x: usize) -> usize {
    // For nonzero x, next_power_of_two() does the job.
    if x == 0 {
        1
    } else {
        x.next_power_of_two()
    }
}

// The stream wrapper. F is the underlying file/stream type.
pub struct KStream<F> {
    pub buf: Vec<u8>,
    pub begin: usize,
    pub end: usize,
    pub is_eof: bool,
    pub f: F,
}

impl<F> KStream<F> {
    /// Initialize a new stream with the given file/stream and buffer size.
    pub fn new(f: F, bufsize: usize) -> Self {
        Self {
            buf: vec![0; bufsize],
            begin: 0,
            end: 0,
            is_eof: false,
            f,
        }
    }

    /// Check if we have reached EOF.
    pub fn eof(&self) -> bool {
        self.is_eof && self.begin >= self.end
    }

    /// Reset the stream.
    pub fn rewind(&mut self) {
        self.is_eof = false;
        self.begin = 0;
        self.end = 0;
    }

    /// Get one byte from the stream.
    pub fn getc(&mut self, mut read: impl FnMut(&mut F, &mut [u8]) -> usize) -> i32 {
        if self.eof() {
            return -1;
        }
        if self.begin >= self.end {
            self.begin = 0;
            self.end = read(&mut self.f, &mut self.buf[..]);
            if self.end == 0 {
                self.is_eof = true;
                return -1;
            }
        }
        let c = self.buf[self.begin];
        self.begin += 1;
        c as i32
    }

    /// Read until a delimiter is encountered. The 'append' flag
    /// indicates whether to clear the target string first.
    pub fn getuntil2(
        &mut self,
        delimiter: i32,
        target: &mut String,
        dret: Option<&mut i32>,
        append: bool,
        mut read: impl FnMut(&mut F, &mut [u8]) -> usize,
    ) -> i32 {
        let mut got_any = false;
        if let Some(dret_ref) = dret {
            *dret_ref = 0;
        }
        if !append {
            target.clear();
        }
        loop {
            if self.begin >= self.end {
                if !self.is_eof {
                    self.begin = 0;
                    self.end = read(&mut self.f, &mut self.buf[..]);
                    if self.end == 0 {
                        self.is_eof = true;
                        break;
                    }
                } else {
                    break;
                }
            }
            let slice = &self.buf[self.begin..self.end];
            let pos = if delimiter == KS_SEP_LINE {
                slice.iter().position(|&b| b == b'\n')
            } else if delimiter > KS_SEP_MAX {
                let d = delimiter as u8;
                slice.iter().position(|&b| b == d)
            } else if delimiter == KS_SEP_SPACE {
                slice.iter().position(|&b| b.is_ascii_whitespace())
            } else if delimiter == KS_SEP_TAB {
                slice
                    .iter()
                    .position(|&b| b.is_ascii_whitespace() && b != b' ')
            } else {
                None
            };
            let i = match pos {
                Some(p) => self.begin + p,
                None => self.end,
            };

            let segment = &self.buf[self.begin..i];
            // Assume valid UTF-8 (or fallback via lossless conversion)
            if let Ok(s) = std::str::from_utf8(segment) {
                target.push_str(s);
            } else {
                target.push_str(&String::from_utf8_lossy(segment));
            }
            got_any = true;
            self.begin = i + 1;
            if i < self.end {
                if let Some(dret_ref) = dret {
                    *dret_ref = self.buf[i] as i32;
                }
                break;
            }
        }
        if !got_any && self.eof() {
            return -1;
        }
        // If reading a line, remove a trailing carriage return.
        if delimiter == KS_SEP_LINE && target.ends_with('\r') {
            target.pop();
        }
        target.len() as i32
    }

    /// Convenience wrapper for getuntil2 without appending.
    pub fn getuntil(
        &mut self,
        delimiter: i32,
        target: &mut String,
        dret: Option<&mut i32>,
        read: impl FnMut(&mut F, &mut [u8]) -> usize,
    ) -> i32 {
        self.getuntil2(delimiter, target, dret, false, read)
    }
}

// Custom error type for sequence reading.
#[derive(Debug)]
pub enum KSeqError {
    Eof,
    TruncatedQuality,
}

// KSeq represents a parsed sequence record (FASTA/FASTQ).
pub struct KSeq<F> {
    pub name: String,
    pub comment: String,
    pub seq: String,
    pub qual: String,
    pub last_char: i32,
    pub f: KStream<F>,
}

impl<F> KSeq<F> {
    /// Initialize a new KSeq using the given file/stream.
    pub fn new(fd: F) -> Self {
        // Use a 16K buffer as in the original.
        let bufsize = 16384;
        Self {
            name: String::new(),
            comment: String::new(),
            seq: String::new(),
            qual: String::new(),
            last_char: 0,
            f: KStream::new(fd, bufsize),
        }
    }

    /// Reset the record parser.
    pub fn rewind(&mut self) {
        self.last_char = 0;
        self.f.rewind();
    }

    /// Read the next sequence record.
    ///
    /// On success, returns Ok(sequence length). On EOF or errors,
    /// returns an appropriate KSeqError.
    pub fn read_record(
        &mut self,
        mut read: impl FnMut(&mut F, &mut [u8]) -> usize,
    ) -> Result<usize, KSeqError> {
        let ks = &mut self.f;
        let mut c: i32;
        if self.last_char == 0 {
            loop {
                c = ks.getc(&mut read);
                if c == -1 || c == ('>' as i32) || c == ('@' as i32) {
                    break;
                }
            }
            if c == -1 {
                return Err(KSeqError::Eof);
            }
            self.last_char = c;
        }
        self.comment.clear();
        self.seq.clear();
        self.qual.clear();

        // Read header line into 'name'
        if ks.getuntil(0, &mut self.name, Some(&mut c), &mut read) < 0 {
            return Err(KSeqError::Eof);
        }
        if c != ('\n' as i32) {
            ks.getuntil(KS_SEP_LINE, &mut self.comment, None, &mut read);
        }
        // Read sequence lines.
        loop {
            c = ks.getc(&mut read);
            if c == -1 || c == ('>' as i32) || c == ('+' as i32) || c == ('@' as i32) {
                break;
            }
            if c == ('\n' as i32) {
                continue;
            }
            self.seq.push(c as u8 as char);
            // Append the rest of the line.
            ks.getuntil2(KS_SEP_LINE, &mut self.seq, None, true, &mut read);
        }
        if c == ('>' as i32) || c == ('@' as i32) {
            self.last_char = c;
        }
        if c != ('+' as i32) {
            // FASTA record.
            return Ok(self.seq.len());
        }
        // FASTQ: read quality header line.
        loop {
            c = ks.getc(&mut read);
            if c == -1 || c == ('\n' as i32) {
                break;
            }
        }
        if c == -1 {
            return Err(KSeqError::TruncatedQuality);
        }
        // Read quality string until it reaches the length of the sequence.
        while ks.getuntil2(KS_SEP_LINE, &mut self.qual, None, true, &mut read) >= 0
            && self.qual.len() < self.seq.len()
        {}
        self.last_char = 0;
        if self.seq.len() != self.qual.len() {
            return Err(KSeqError::TruncatedQuality);
        }
        Ok(self.seq.len())
    }
}
