/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::collections::HashMap;

pub fn translate_to_all_frames(dna_seq: &str, aa_seqs: &mut Vec<String>) {
    // Ensure aa_seqs has 6 elements, one for each reading frame
    aa_seqs.resize(6, String::new());
    let max_size = (dna_seq.len() / 3) + 1;
    for aa_seq in aa_seqs.iter_mut() {
        aa_seq.reserve(max_size);
    }
    if dna_seq.len() < 3 {
        return;
    }

    // Static initialization of translation map and lookup tables
    static INIT: std::sync::Once = std::sync::Once::new();
    static mut FWD_LOOKUP_TABLE: [u8; 256] = [0; 256];
    static mut REV_LOOKUP_TABLE: [u8; 256] = [0; 256];
    static TRANSLATION_MAP: &[u8] =
        b"KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

    INIT.call_once(|| {
        unsafe {
            // Initialize lookup tables with UINT8_MAX (255)
            for i in 0..256 {
                FWD_LOOKUP_TABLE[i] = 255;
                REV_LOOKUP_TABLE[i] = 255;
            }
            // Map is based on AGCT coding, not ACGT
            FWD_LOOKUP_TABLE[b'A' as usize] = 0x00;
            FWD_LOOKUP_TABLE[b'a' as usize] = 0x00;
            FWD_LOOKUP_TABLE[b'G' as usize] = 0x01;
            FWD_LOOKUP_TABLE[b'g' as usize] = 0x01;
            FWD_LOOKUP_TABLE[b'C' as usize] = 0x02;
            FWD_LOOKUP_TABLE[b'c' as usize] = 0x02;
            FWD_LOOKUP_TABLE[b'T' as usize] = 0x03;
            FWD_LOOKUP_TABLE[b't' as usize] = 0x03;

            REV_LOOKUP_TABLE[b'A' as usize] = 0x30;
            REV_LOOKUP_TABLE[b'a' as usize] = 0x30;
            REV_LOOKUP_TABLE[b'G' as usize] = 0x20;
            REV_LOOKUP_TABLE[b'g' as usize] = 0x20;
            REV_LOOKUP_TABLE[b'C' as usize] = 0x10;
            REV_LOOKUP_TABLE[b'c' as usize] = 0x10;
            REV_LOOKUP_TABLE[b'T' as usize] = 0x00;
            REV_LOOKUP_TABLE[b't' as usize] = 0x00;
        }
    });

    let mut fwd_codon: u8 = 0;
    let mut rev_codon: u8 = 0;
    let mut ambig_nt_countdown = 0; // If positive, bases to go until 'N' leaves codon
    let mut frame_len = [0usize; 6];

    let dna_bytes = dna_seq.as_bytes();
    for (i, &nt) in dna_bytes.iter().enumerate() {
        let frame = i % 3;
        unsafe {
            fwd_codon = (fwd_codon << 2) & 0x3F;
            rev_codon >>= 2;

            if ambig_nt_countdown > 0 {
                ambig_nt_countdown -= 1;
            }
            let fwd_lookup_code = FWD_LOOKUP_TABLE[nt as usize];
            let rev_lookup_code = REV_LOOKUP_TABLE[nt as usize];

            if fwd_lookup_code == 255 {
                // Unknown nucleotide, set ambiguity countdown
                ambig_nt_countdown = 3;
            } else {
                fwd_codon |= fwd_lookup_code;
                rev_codon |= rev_lookup_code;
            }

            if i >= 2 {
                // We have a full codon
                let ch_fwd = if ambig_nt_countdown == 0 {
                    TRANSLATION_MAP[fwd_codon as usize] as char
                } else {
                    'X'
                };
                aa_seqs[frame].push(ch_fwd);
                frame_len[frame] += 1;

                let ch_rev = if ambig_nt_countdown == 0 {
                    TRANSLATION_MAP[rev_codon as usize] as char
                } else {
                    'X'
                };
                // Prepend to reverse frame (frames 3 to 5)
                aa_seqs[frame + 3].insert(0, ch_rev);
                frame_len[frame + 3] += 1;
            }
        }
    }

    // Resize aa_seqs to actual lengths
    for i in 0..3 {
        aa_seqs[i].truncate(frame_len[i]);
    }
    for i in 3..6 {
        aa_seqs[i].truncate(frame_len[i]);
    }
}
