use std::sync::OnceLock;

/// A 64-character lookup for codon -> amino acid.
/// Matches the order used in the C++ code: "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF"
/// This covers 4^3 = 64 possible codons in AGCT ordering.
static TRANSLATION_MAP: &str = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

/// Lazy-initialized lookup table for the **forward** strand.
/// Maps ASCII 'A', 'G', 'C', 'T' (uppercase or lowercase) to 0..3 in an AGCT order.
/// All other bytes map to `u8::MAX` to indicate “ambiguous.”
static FWD_LOOKUP_TABLE: OnceLock<[u8; 256]> = OnceLock::new();

/// Lazy-initialized lookup table for the **reverse** strand.
/// Maps 'A', 'G', 'C', 'T' to 0..3 in reverse-complement style, shifted to upper bits.
static REV_LOOKUP_TABLE: OnceLock<[u8; 256]> = OnceLock::new();

/// Initialize the forward and reverse lookup tables on first use.
fn init_lookup_tables() {
    // If already initialized, return.
    if FWD_LOOKUP_TABLE.get().is_some() && REV_LOOKUP_TABLE.get().is_some() {
        return;
    }

    let mut fwd_table = [u8::MAX; 256];
    let mut rev_table = [u8::MAX; 256];

    // 'A' (or 'a') -> 0, 'G' -> 1, 'C' -> 2, 'T' -> 3 in forward
    // Reverse complement is in the top nibble, so we shift accordingly:
    //   'A' -> 0x30, 'G' -> 0x20, 'C' -> 0x10, 'T' -> 0x00 for reverse
    for &base in b"AGCTagct" {
        let code_forward = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'G' => 1,
            b'C' => 2,
            b'T' => 3,
            _ => u8::MAX,
        };
        let code_reverse = match base.to_ascii_uppercase() {
            b'A' => 0x30,
            b'G' => 0x20,
            b'C' => 0x10,
            b'T' => 0x00,
            _ => u8::MAX,
        };
        fwd_table[base as usize] = code_forward;
        rev_table[base as usize] = code_reverse;
    }

    let _ = FWD_LOOKUP_TABLE.set(fwd_table);
    let _ = REV_LOOKUP_TABLE.set(rev_table);
}

/// Translate a single codon (0..63) using `TRANSLATION_MAP` or returns `'X'` if ambiguous.
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
pub fn translate_to_all_frames(dna_seq: &str) -> [String; 6] {
    init_lookup_tables();
    let fwd_table = FWD_LOOKUP_TABLE.get().unwrap();
    let rev_table = REV_LOOKUP_TABLE.get().unwrap();

    // We allocate 6 frames; they might not all be full, so we'll collect them in vectors of chars.
    // We'll finalize them as Strings at the end.
    let mut aa_frames = [
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ];

    // If the sequence is shorter than 3, all frames remain empty.
    if dna_seq.len() < 3 {
        return [
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        ];
    }

    // We'll track the forward codon bits (0..63) and reverse codon bits (0..63).
    let mut fwd_codon: u8 = 0;
    let mut rev_codon: u8 = 0;

    // If ambig_nt_countdown > 0, we are in a window of 3 bases that has at least one ambiguous char.
    let mut ambig_nt_countdown: u8 = 0;

    for (i, &nt) in dna_seq.as_bytes().iter().enumerate() {
        let frame = i % 3;

        // Shift forward codon left by 2 bits, keeping only lower 6 bits
        fwd_codon = (fwd_codon << 2) & 0x3F;
        // Shift reverse codon right by 2 bits, discarding lower bits
        rev_codon >>= 2;

        // Decrease ambiguous countdown if active
        if ambig_nt_countdown > 0 {
            ambig_nt_countdown -= 1;
        }

        // Lookup forward and reverse codes
        let fwd_lookup_code = fwd_table[nt as usize];
        let rev_lookup_code = rev_table[nt as usize];

        if fwd_lookup_code == u8::MAX {
            // Mark 3 bases as ambiguous
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_lookup_code;
            rev_codon |= rev_lookup_code;
        }

        // Once we have read at least 3 bases, we have a full codon
        if i >= 2 {
            // Translate forward codon
            let ch_fwd = translate_codon(fwd_codon, ambig_nt_countdown > 0);
            aa_frames[frame].push(ch_fwd);

            // Translate reverse codon
            let ch_rev = translate_codon(rev_codon, ambig_nt_countdown > 0);
            // The original code prepends the amino acid (building from the end),
            // but we can build from the front and reverse later
            aa_frames[frame + 3].push(ch_rev);
        }
    }

    // Now, we need to reverse the 3 reverse-frame strings because in the C++ code,
    // they were filled from the “end” backwards.
    for frame in 3..6 {
        aa_frames[frame].reverse();
    }

    // Convert from Vec<char> to String
    let mut result = [
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
    ];
    for (i, frame_vec) in aa_frames.into_iter().enumerate() {
        result[i] = frame_vec.into_iter().collect();
    }

    result
}
