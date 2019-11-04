#[cfg(test)]
mod tests {
    use super::hamming::*;
    use super::*;
    #[test]
    fn test_hamming_str_dist() {
        assert_eq!(hamming_str("ACTG", "ACTT"), 1);
        assert_eq!(hamming_str("ACTG", "ACTTT"), 1);
    }

    #[test]
    fn test_base_bits() {
        let alpha = BaseBits::new(b"ACTG").unwrap();
        let beta = BaseBits::new(b"ACTT").unwrap();
        assert_eq!(hamming_dist(&alpha, &beta), 1);
    }

    #[test]
    fn test_cases_bb_hamming() {
        // Test N encoding
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"ACTN").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );
        // Test * encoding
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"ACT*").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );
        // Test that UNDETERMINED encoding
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"ACT9").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            1
        );
        // Test regular equality
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"ACTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );

        // Test Other string
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"GATACA").unwrap(),
                &BaseBits::new(b"GATACT").unwrap()
            ),
            1
        );
        assert_eq!(
            hamming_dist(
                &BaseBits::new(b"CATACAGATACTTCCATAGCT").unwrap(),
                &BaseBits::new(b"GATACAGATACTTCCATAGCA").unwrap()
            ),
            2
        );
        // This one should fail since it overflows and wraps around. needs thought
        //assert_eq!(hamming_dist(&BaseBits::new("TATACAGATACTTCCATAGCATC"),
        //&BaseBits::new("GATACAGATACAACNATAGCATT")), 4);
    }

    #[test]
    fn test_bb_to_string() {
        let alpha = BaseBits::new(b"GCTAN").unwrap();
        let beta = BaseBits::new(b"ACTG*").unwrap();
        assert_eq!(alpha.to_string(), "GCTAN".to_string());
        assert_eq!(beta.to_string(), "ACTGN".to_string());

        let long = BaseBits::new(b"GATACAGATACAACNATAGCA").unwrap();
        assert_eq!(long.to_string(), "GATACAGATACAACNATAGCA".to_string());
    }

    #[test]
    fn test_encoding() {
        let bb = BaseBits::new(b"ACTG").unwrap();
        assert_eq!(bb.code, 0b000110101011);
    }
}

/// Encode a DNA string of up to 21 bases as a u64 for fast hamming distance calculations.
/// TODO: Add a bump to use u128 or maybe bigint if 21 chars is not enough.
/// TODO: Add equalities and hash function stuff so this type can be used in data structures
use std::fmt;
use std::str;
use std::u64;

pub const ENCODING_DIST: u32 = 2;
pub const ENCODING_LENGTH: u32 = 3;
pub const CONTAINER_WIDTH: u32 = 64;
pub const MAX_BASES: usize = (CONTAINER_WIDTH / ENCODING_LENGTH) as usize;
pub const UNDETERMINED: u64 = 0b100;
pub const ANY: u64 = 0b111;
pub const MAX_VAL: u64 = u64::MAX;

struct Bases;
impl Bases {
    const A: u64 = 0b000;
    const C: u64 = 0b110;
    const T: u64 = 0b101;
    const G: u64 = 0b011;
    const N: u64 = ANY;
    const STAR: u64 = ANY;
    const UNDETERMINED: u64 = UNDETERMINED;
}

//#[derive(Copy, Clone)]
#[derive(Hash, PartialEq, Eq, Debug, Copy, Clone)]
pub struct BaseBits {
    pub code: u64,
    len: usize,
}

impl BaseBits {
    pub fn new(seq: &[u8]) -> Result<BaseBits, &'static str> {
        let mut code: u64 = 0;
        let len = seq.len();
        if len > MAX_BASES {
            return Err("Length of string to encode exceeds MAX_BASES");
        }
        for c in seq.iter() {
            code = (code << ENCODING_LENGTH)
                | match c {
                    b'A' => Bases::A,
                    b'C' => Bases::C,
                    b'T' => Bases::T,
                    b'G' => Bases::G,
                    b'N' => Bases::N,
                    b'*' => Bases::STAR,
                    _ => Bases::UNDETERMINED,
                }
        }
        Ok(BaseBits { code, len })
    }

    //#[allow(unreachable_patterns)]
    pub fn decode(&self) -> Vec<u8> {
        let mut s = Vec::new();
        let mut code = self.code;
        for _ in 0..self.len {
            let base = extract_bits(code, ENCODING_LENGTH);
            code = code >> ENCODING_LENGTH;
            s.push(match base {
                Bases::A => b'A',
                Bases::C => b'C',
                Bases::T => b'T',
                Bases::G => b'G',
                _ => b'N',
                //Bases::N => b'N',
                //Bases::STAR => b'*', // As long as N and STAR are ANY, STAR is unreachable
                //_ => b'N',
            });
        }
        s.into_iter().rev().collect()
    }
}

impl fmt::Display for BaseBits {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", str::from_utf8(&self.decode()).unwrap())
    }
}

#[inline]
pub fn hamming_dist(alpha: &BaseBits, beta: &BaseBits) -> u32 {
    (alpha.code ^ beta.code).count_ones() / 2
}

/// Extract 'k' bits from the end of a u64 integer
#[inline]
fn extract_bits(n: u64, k: u32) -> u64 {
    !(!0u64 << k) & n
}

// Hamming distance functions that don't depend on BaseBits types
pub mod hamming {
    #[inline]
    pub fn hamming_code(alpha: u64, beta: u64) -> u32 {
        //let x = alpha ^ beta;
        //x.count_ones() / 2
        (alpha ^ beta).count_ones() / 2
    }

    pub fn hamming_str(alpha: &str, beta: &str) -> u32 {
        // skip length check. will default to up to length of alpha
        let mut dist = 0;
        for (a, b) in alpha.chars().zip(beta.chars()) {
            if a != b {
                dist += 1
            }
        }
        dist
    }
}

// dist (pop_count in the paper)
// https://github.com/Daniel-Liu-c0deb0t/UMICollapse/blob/master/src/umicollapse/util/Read.java
