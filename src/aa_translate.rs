use std::sync::OnceLock;

/// A 64-character lookup for codon -> amino acid.
/// Matches the order used in the C++ code: "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF"
/// This covers 4^3 = 64 possible codons in AGCT ordering.
static TRANSLATION_MAP: &str = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

/// Lazy-initialized lookup table for the **forward** strand.
/// Maps ASCII 'A', 'G', 'C', 'T' (uppercase or lowercase) to 0..3 in an AGCT order.
/// All other bytes map to `u8::MAX` to indicate "ambiguous."
static FWD_LOOKUP_TABLE: OnceLock<[u8; 256]> = OnceLock::new();

/// Initialize the forward lookup table on first use.
fn init_lookup_table() {
    // If already initialized, return.
    if FWD_LOOKUP_TABLE.get().is_some() {
        return;
    }

    let mut fwd_table = [u8::MAX; 256];

    // 'A' (or 'a') -> 0, 'G' -> 1, 'C' -> 2, 'T' -> 3 in forward
    for &base in b"AGCTagct" {
        let code = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'G' => 1,
            b'C' => 2,
            b'T' => 3,
            _ => u8::MAX,
        };
        fwd_table[base as usize] = code;
    }

    let _ = FWD_LOOKUP_TABLE.set(fwd_table);
}

/// Translate a single codon (packed into a u8) into an amino acid char.
#[inline]
fn translate_codon(codon: u8, is_ambiguous: bool) -> char {
    if is_ambiguous {
        'X'
    } else {
        // Index into the translation map.
        let idx = codon as usize & 0x3F; // ensure it's 0..63
        TRANSLATION_MAP.as_bytes()[idx] as char
    }
}

/// Translates a DNA sequence into all 6 reading frames (3 forward, 3 reverse).
///
/// The returned array `[String; 6]` contains:
/// - Index 0..2: forward frames
/// - Index 3..5: reverse frames
///
/// # Parameters
/// - `dna_seq`: The DNA sequence. Must be uppercase/lowercase letters A/C/G/T (others treated as ambiguous).
///
/// # Returns
/// An array of 6 amino-acid sequences (forward frames 0..2, reverse frames 3..5).
///
pub fn translate_to_all_frames(seq: &str) -> Vec<String> {
    let mut frames = vec![String::new(); 6];
    // Forward frames
    for frame in 0..3 {
        translate_frame(seq, frame, &mut frames[frame]);
    }

    // Reverse frames
    let rev_seq = reverse_complement(seq);
    for frame in 0..3 {
        translate_frame(&rev_seq, frame, &mut frames[frame + 3]);
    }
    frames
}

fn translate_frame(seq: &str, frame: usize, output: &mut String) {
    init_lookup_table(); // Make sure lookup table is initialized
    let bases = seq.as_bytes();
    let mut i = frame;
    while i + 2 < bases.len() {
        let codon = &bases[i..i + 3];
        let is_ambiguous = codon
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));
        let packed = pack_codon(codon);
        output.push(translate_codon(packed, is_ambiguous));
        i += 3;
    }
}

fn pack_codon(codon: &[u8]) -> u8 {
    let mut packed = 0u8;
    let lookup = FWD_LOOKUP_TABLE.get().unwrap();
    for (i, &base) in codon.iter().enumerate() {
        let code = lookup[base as usize];
        packed |= code << (2 * (2 - i));
    }
    packed
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            _ => 'N',
        })
        .collect()
}
