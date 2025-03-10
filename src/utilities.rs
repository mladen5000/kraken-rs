// Define CURRENT_REVCOM_VERSION from build_db
pub const CURRENT_REVCOM_VERSION: i32 = 2;

pub fn expand_spaced_seed_mask(mask: &mut u64, bits_per_char: i32) {
    let mut new_mask = 0u64;
    let mut pos = 0;
    let mut old_pos = 0;

    while old_pos < 64 {
        if (*mask & (1u64 << old_pos)) != 0 {
            new_mask |= 1u64 << pos;
            pos += bits_per_char as u64;
        }
        old_pos += 1;
    }
    *mask = new_mask;
}

pub fn split_string(str_val: &str, delim: &str, max_fields: usize) -> Vec<String> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut field_count = 0;

    while field_count < max_fields - 1 {
        if let Some(pos) = str_val[start..].find(delim) {
            fields.push(str_val[start..start + pos].to_string());
            start += pos + delim.len();
            field_count += 1;
        } else {
            break;
        }
    }

    if start < str_val.len() {
        fields.push(str_val[start..].to_string());
    }

    fields
}

pub fn murmur_hash3(k: u64) -> u64 {
    let mut h = k;

    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;

    h
}

pub fn murmur_hash3_64(bytes: &[u8]) -> u64 {
    let len = bytes.len();
    let seed = 0;
    let mut h1: u64 = seed;
    let mut h2: u64 = seed;

    let c1: u64 = 0x87c37b91114253d5;
    let c2: u64 = 0x4cf5ad432745937f;

    // Body
    let nblocks = len / 16;
    let blocks = unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, nblocks * 2) };

    for i in 0..nblocks {
        let mut k1 = blocks[i * 2];
        let mut k2 = blocks[i * 2 + 1];

        k1 = k1.wrapping_mul(c1);
        k1 = k1.rotate_left(31);
        k1 = k1.wrapping_mul(c2);
        h1 ^= k1;

        h1 = h1.rotate_left(27);
        h1 = h1.wrapping_add(h2);
        h1 = h1.wrapping_mul(5).wrapping_add(0x52dce729);

        k2 = k2.wrapping_mul(c2);
        k2 = k2.rotate_left(33);
        k2 = k2.wrapping_mul(c1);
        h2 ^= k2;

        h2 = h2.rotate_left(31);
        h2 = h2.wrapping_add(h1);
        h2 = h2.wrapping_mul(5).wrapping_add(0x38495ab5);
    }

    // Tail
    let tail = &bytes[nblocks * 16..];
    let mut k1: u64 = 0;
    let mut k2: u64 = 0;

    match tail.len() {
        15 => {
            k2 ^= (tail[14] as u64) << 48;
            k2 ^= (tail[13] as u64) << 40;
            k2 ^= (tail[12] as u64) << 32;
            k2 ^= (tail[11] as u64) << 24;
            k2 ^= (tail[10] as u64) << 16;
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        14 => {
            k2 ^= (tail[13] as u64) << 40;
            k2 ^= (tail[12] as u64) << 32;
            k2 ^= (tail[11] as u64) << 24;
            k2 ^= (tail[10] as u64) << 16;
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        13 => {
            k2 ^= (tail[12] as u64) << 32;
            k2 ^= (tail[11] as u64) << 24;
            k2 ^= (tail[10] as u64) << 16;
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        12 => {
            k2 ^= (tail[11] as u64) << 24;
            k2 ^= (tail[10] as u64) << 16;
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        11 => {
            k2 ^= (tail[10] as u64) << 16;
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        10 => {
            k2 ^= (tail[9] as u64) << 8;
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        9 => {
            k2 ^= tail[8] as u64;
            k2 = k2.wrapping_mul(c2);
            k2 = k2.rotate_left(33);
            k2 = k2.wrapping_mul(c1);
            h2 ^= k2;

            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        8 => {
            k1 ^= (tail[7] as u64) << 56;
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        7 => {
            k1 ^= (tail[6] as u64) << 48;
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        6 => {
            k1 ^= (tail[5] as u64) << 40;
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        5 => {
            k1 ^= (tail[4] as u64) << 32;
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        4 => {
            k1 ^= (tail[3] as u64) << 24;
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        3 => {
            k1 ^= (tail[2] as u64) << 16;
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        2 => {
            k1 ^= (tail[1] as u64) << 8;
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        1 => {
            k1 ^= tail[0] as u64;
            k1 = k1.wrapping_mul(c1);
            k1 = k1.rotate_left(31);
            k1 = k1.wrapping_mul(c2);
            h1 ^= k1;
        }
        _ => {}
    }

    // Finalization
    h1 ^= len as u64;
    h2 ^= len as u64;

    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 = h1.wrapping_add(h2);

    h1
}

fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_spaced_seed_mask() {
        let mut mask = 0b010110u64;
        expand_spaced_seed_mask(&mut mask, 2);
        // Each 1 bit in input expands to consecutive 1s in output
        // Each 0 bit in input expands to consecutive 0s in output
        // The number of consecutive bits is determined by bits_per_char
        let expected = 0b001100001100110000u64;
        assert_eq!(mask, expected);
    }

    #[test]
    fn test_split_string() {
        let s = "field1\tfield2\tfield3\tfield4";
        let result = split_string(s, "\t", 3);
        assert_eq!(result, vec!["field1", "field2", "field3\tfield4"]);
    }
}
