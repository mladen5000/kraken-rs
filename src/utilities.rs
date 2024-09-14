/// Expands a spaced seed mask based on a bit expansion factor.
///
/// # Arguments
///
/// * `spaced_seed_mask` - The original spaced seed mask to be expanded.
/// * `bit_expansion_factor` - The factor by which each bit in the mask is expanded.
///
/// # Returns
///
/// The expanded spaced seed mask.
pub fn expand_spaced_seed_mask(spaced_seed_mask: &mut u64, bit_expansion_factor: usize) {
    let mut new_mask: u64 = 0;
    let bits = (1 << bit_expansion_factor) - 1;

    for i in 0..64 {
        new_mask <<= bit_expansion_factor;
        if (*spaced_seed_mask >> (63 - i)) & 1 != 0 {
            new_mask |= bits;
        }
    }
    *spaced_seed_mask = new_mask;
}

/// Splits a string based on a delimiter into a vector of substrings.
///
/// # Arguments
///
/// * `str` - The string to be split.
/// * `delim` - The delimiter to split the string.
/// * `max_fields` - The maximum number of fields to split into.
///
/// # Returns
///
/// A vector containing the split substrings.
pub fn split_string(str: &str, delim: &str, max_fields: usize) -> Vec<String> {
    let mut output = Vec::new();
    let mut pos1 = 0;
    let mut field_ct = 0;
    let mut finished = false;

    while field_ct < max_fields - 1 && !finished {
        let pos2 = str[pos1..].find(delim).map(|x| x + pos1);
        let token = match pos2 {
            None => {
                finished = true;
                str[pos1..].to_string()
            }
            Some(pos) => {
                let token = &str[pos1..pos];
                pos1 = pos + delim.len();
                token.to_string()
            }
        };
        output.push(token);
        field_ct += 1;
    }

    if pos1 < str.len() {
        output.push(str[pos1..].to_string());
    }

    output
}

fn main() {
    // Example usage of expand_spaced_seed_mask
    let mut mask: u64 = 0b1010;
    expand_spaced_seed_mask(&mut mask, 3);
    println!("Expanded mask: {:064b}", mask);

    // Example usage of split_string
    let input_str = "one,two,three,four";
    let split = split_string(input_str, ",", 3);
    println!("Split string: {:?}", split);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_spaced_seed_mask() {
        let mut mask: u64 = 0b1010;
        expand_spaced_seed_mask(&mut mask, 3);
        assert_eq!(mask, 0b111000000111000000);
    }

    #[test]
    fn test_split_string() {
        let input_str = "one,two,three,four";
        let split = split_string(input_str, ",", 3);
        assert_eq!(split, vec!["one", "two", "three,four"]);

        let split_all = split_string(input_str, ",", 10);
        assert_eq!(split_all, vec!["one", "two", "three", "four"]);
    }
}
