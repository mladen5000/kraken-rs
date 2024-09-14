/*
 * Copyright 2013-2023, Derrick Wood
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

use std::u64;

/// Expands a spaced seed mask by repeating each bit `bit_expansion_factor` times.
///
/// # Arguments
///
/// * `spaced_seed_mask` - A mutable reference to the original seed mask.
/// * `bit_expansion_factor` - The number of times to repeat each bit.
///
/// # Example
///
/// ```
/// let mut mask = 0b101;
/// expand_spaced_seed_mask(&mut mask, 2);
/// assert_eq!(mask, 0b11_00_11);
/// ```
pub fn expand_spaced_seed_mask(spaced_seed_mask: &mut u64, bit_expansion_factor: i32) {
    let bits = (1u64 << bit_expansion_factor) - 1;
    let mut new_mask: u64 = 0;

    let n = 64 / bit_expansion_factor;

    for i in (0..n).rev() {
        new_mask <<= bit_expansion_factor as u32;
        if ((*spaced_seed_mask >> i as u32) & 1) != 0 {
            new_mask |= bits;
        }
    }
    *spaced_seed_mask = new_mask;
}

/// Splits a string into a vector of substrings based on a delimiter, up to a maximum number of fields.
///
/// # Arguments
///
/// * `s` - The string to split.
/// * `delim` - The delimiter string.
/// * `max_fields` - The maximum number of fields to split into.
///
/// # Returns
///
/// A vector of substrings.
///
/// # Example
///
/// ```
/// let s = "one,two,three,four";
/// let result = split_string(&s, ",", 3);
/// assert_eq!(result, vec!["one", "two", "three,four"]);
/// ```
pub fn split_string(s: &str, delim: &str, max_fields: usize) -> Vec<String> {
    let mut output = Vec::new();
    let mut pos1 = 0;
    let mut field_ct = 0;
    let mut finished = false;

    while field_ct < max_fields && !finished {
        if let Some(pos2_rel) = s[pos1..].find(delim) {
            let pos2 = pos1 + pos2_rel;
            let token = &s[pos1..pos2];
            output.push(token.to_string());
            pos1 = pos2 + delim.len();
        } else {
            let token = &s[pos1..];
            output.push(token.to_string());
            finished = true;
        }
        field_ct += 1;
    }
    output
}
