use std::sync::Once;

// Translation map for codons to amino acids.
const TRANSLATION_MAP: &[char] = &[
    'K', 'K', 'N', 'N', 'R', 'R', 'S', 'S', 'T', 'T', 'T', 'T', 'I', 'M', 'I', 'I', 'E', 'E', 'D',
    'D', 'G', 'G', 'G', 'G', 'A', 'A', 'A', 'A', 'V', 'V', 'V', 'V', 'Q', 'Q', 'H', 'H', 'R', 'R',
    'R', 'R', 'P', 'P', 'P', 'P', 'L', 'L', 'L', 'L', '*', '*', 'Y', 'Y', '*', 'W', 'C', 'C', 'S',
    'S', 'S', 'S', 'L', 'L', 'F', 'F',
];

// Lookup tables for DNA nucleotides.
static INIT: Once = Once::new();
static mut FWD_LOOKUP_TABLE: [Option<u8>; 256] = [None; 256];
static mut REV_LOOKUP_TABLE: [Option<u8>; 256] = [None; 256];

/// Builds the forward and reverse lookup tables for DNA nucleotides.
fn build_lookup_tables() {
    unsafe {
        for &ch in &['A', 'G', 'C', 'T', 'a', 'g', 'c', 't'] {
            let index = match ch {
                'A' | 'a' => Some(0x00),
                'G' | 'g' => Some(0x01),
                'C' | 'c' => Some(0x02),
                'T' | 't' => Some(0x03),
                _ => unreachable!(),
            };
            FWD_LOOKUP_TABLE[ch as usize] = index;
            // The reverse lookup table is built by flipping the bits of the forward table.
            REV_LOOKUP_TABLE[ch as usize] = index.map(|x| x << 4);
        }
    }
}

/// Translates a DNA sequence to amino acid sequences in all six reading frames.
///
/// # Arguments
///
/// * `dna_seq` - The DNA sequence to be translated.
///
/// # Returns
///
/// A vector containing amino acid sequences for each of the six reading frames.
pub fn translate_to_all_frames(dna_seq: &str) -> Vec<String> {
    let max_size = (dna_seq.len() / 3) + 1;
    let mut aa_seqs = vec![String::with_capacity(max_size); 6];
    if dna_seq.len() < 3 {
        return aa_seqs;
    }

    INIT.call_once(|| build_lookup_tables());

    let mut fwd_codon = 0u8;
    let mut rev_codon = 0u8;
    let mut ambig_nt_countdown = 0;
    let mut frame_len = [0usize; 6];

    for (i, nt) in dna_seq.bytes().enumerate() {
        let frame = i % 3;
        fwd_codon = (fwd_codon << 2) & 0x3f;
        rev_codon >>= 2;

        let fwd_lookup_code;
        let rev_lookup_code;
        unsafe {
            fwd_lookup_code = FWD_LOOKUP_TABLE[nt as usize];
            rev_lookup_code = REV_LOOKUP_TABLE[nt as usize];
        }

        if fwd_lookup_code.is_none() {
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_lookup_code.unwrap();
            rev_codon |= rev_lookup_code.unwrap();
        }

        if i >= 2 && ambig_nt_countdown == 0 {
            // Translate and append to forward frame
            aa_seqs[frame].push(TRANSLATION_MAP[fwd_codon as usize]);
            frame_len[frame] += 1;

            // Translate and append to reverse frame
            aa_seqs[frame + 3].insert(0, TRANSLATION_MAP[rev_codon as usize]);
            frame_len[frame + 3] += 1;
        }
    }

    // Truncate sequences to actual lengths
    for i in 0..6 {
        aa_seqs[i].truncate(frame_len[i]);
    }

    aa_seqs
}

fn main() {
    let dna_seq = "AGCTAGCTAGCT";
    let translated_frames = translate_to_all_frames(dna_seq);
    for (i, frame) in translated_frames.iter().enumerate() {
        println!("Frame {}: {}", i + 1, frame);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_to_all_frames() {
        let dna_seq = "AGCTAGCTAGCT";
        let translated_frames = translate_to_all_frames(dna_seq);
        println!("{translated_frames:?}");
        assert_eq!(translated_frames.len(), 6);
        // Example check for specific frame translation (can be expanded for more detailed tests)
        assert_eq!(translated_frames[0], "SSS");
        assert_eq!(translated_frames[5], "LCL");
    }
}
