const TRANSLATION_MAP: [u8; 64] = [
    b'K', b'K', b'N', b'N', b'R', b'R', b'S', b'S', b'T', b'T', b'T', b'T', b'I', b'M', b'I', b'I',
    b'E', b'E', b'D', b'D', b'G', b'G', b'G', b'G', b'A', b'A', b'A', b'A', b'V', b'V', b'V', b'V',
    b'Q', b'Q', b'H', b'H', b'R', b'R', b'R', b'R', b'P', b'P', b'P', b'P', b'L', b'L', b'L', b'L',
    b'*', b'*', b'Y', b'Y', b'*', b'W', b'C', b'C', b'S', b'S', b'S', b'S', b'L', b'L', b'F', b'F',
];

const FWD_LOOKUP_TABLE: [u8; 256] = {
    let mut table = [u8::MAX; 256];
    table[b'A' as usize] = 0x00;
    table[b'G' as usize] = 0x01;
    table[b'C' as usize] = 0x02;
    table[b'T' as usize] = 0x03;
    table
};

const REV_LOOKUP_TABLE: [u8; 256] = {
    let mut table = [u8::MAX; 256];
    table[b'A' as usize] = 0x30;
    table[b'G' as usize] = 0x20;
    table[b'C' as usize] = 0x10;
    table[b'T' as usize] = 0x00;
    table
};

fn translate_to_all_frames(dna_seq: &str) -> [String; 6] {
    let max_size = (dna_seq.len() / 3) + 1;
    let mut aa_seqs = [
        Vec::with_capacity(max_size),
        Vec::with_capacity(max_size),
        Vec::with_capacity(max_size),
        Vec::with_capacity(max_size),
        Vec::with_capacity(max_size),
        Vec::with_capacity(max_size),
    ];

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

    let mut fwd_codon = 0u8;
    let mut rev_codon = 0u8;
    let mut ambig_nt_countdown = 0;

    for (i, base) in dna_seq.bytes().enumerate() {
        let frame = i % 3;
        fwd_codon = (fwd_codon << 2) & 0x3f;
        rev_codon >>= 2;
        if ambig_nt_countdown > 0 {
            ambig_nt_countdown -= 1;
        }

        let fwd_lookup_code = FWD_LOOKUP_TABLE[base as usize];
        let rev_lookup_code = REV_LOOKUP_TABLE[base as usize];
        if fwd_lookup_code == u8::MAX {
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_lookup_code;
            rev_codon |= rev_lookup_code;
        }

        if i >= 2 {
            let fwd_ch = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[fwd_codon as usize]
            } else {
                b'X'
            };
            aa_seqs[frame].push(fwd_ch);

            let rev_ch = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[rev_codon as usize]
            } else {
                b'X'
            };
            aa_seqs[frame + 3].push(rev_ch);
        }
    }

    for seq in &mut aa_seqs[3..6] {
        seq.reverse();
    }

    let aa_seqs: [String; 6] = aa_seqs
        .into_iter()
        .map(|vec| unsafe { String::from_utf8_unchecked(vec) })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    aa_seqs
}
