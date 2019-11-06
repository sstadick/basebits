/// Encode a DNA string of up to 21 bases as a u64 for fast hamming distance calculations.
/// Each BaseBits will take up u64 X 2 + usize amount of space. It works by having encodings for A,
/// C, T, and G that are all dist 2 away from eachother. A sequence is encoded into a u64 by
/// setting the bits for each character. Any unrecognized character is treated as an N. N's are
/// encoded as 001, but are also tracked speratalty to allow for two different methods of counting.
/// N's can be treated as wildcards by using the `hamming_dist_nany` method, or they can be treated
/// like a character by using the `hamming_dist_none` method.
///
/// Generally speaking, if you are going to compare against a string more than 4 times, it is worth
/// the cost of encoding it and using this package.
use std::fmt;
use std::str;
use std::u64;

pub const ENCODING_DIST: u32 = 2;
pub const ENCODING_LENGTH: u32 = 3;
pub const CONTAINER_WIDTH: u32 = 64;
pub const MAX_BASES: usize = (CONTAINER_WIDTH / ENCODING_LENGTH) as usize;
pub const UNDETERMINED: u64 = 0b100;
//pub const ANY: u64 = 0b111;
pub const MAX_VAL: u64 = u64::MAX;

struct Bases;
impl Bases {
    const A: u64 = 0b000;
    const C: u64 = 0b110;
    const T: u64 = 0b101;
    const G: u64 = 0b011;
    const N: u64 = UNDETERMINED;
}

/// A BaseBits encoding
#[derive(Hash, PartialEq, Eq, Debug, Copy, Clone)]
pub struct BaseBits {
    /// The u64 holding the encoding
    pub code: u64,
    /// The u64 holding an inverse encoding of N's
    nbits: u64,
    /// The length of the original input
    len: usize,
}

impl BaseBits {
    /// Create a new BaseBits object.
    pub fn new(seq: &[u8]) -> Result<BaseBits, &'static str> {
        let mut code: u64 = 0;
        let mut nbits: u64 = !0b0;
        let len = seq.len();
        if len > MAX_BASES {
            return Err("Length of string to encode exceeds MAX_BASES");
        }
        for c in seq.iter() {
            let base = match c {
                b'A' => Bases::A,
                b'C' => Bases::C,
                b'T' => Bases::T,
                b'G' => Bases::G,
                _ => Bases::N,
            };

            code = (code << ENCODING_LENGTH) | base;
            nbits = match base {
                Bases::N => (nbits << ENCODING_LENGTH) | 0b000,
                _ => (nbits << ENCODING_LENGTH) | 0b111,
            }
        }
        Ok(BaseBits { code, nbits, len })
    }

    /// Decode a BaseBits object into a string
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

/// Compute hamming distance between two strings, count N's as any character
#[inline]
pub fn hamming_dist_nany(alpha: &BaseBits, beta: &BaseBits) -> u32 {
    ((alpha.code ^ beta.code) & (alpha.nbits & beta.nbits)).count_ones() / ENCODING_DIST
}

/// Compute hamming distace but N's as mismatches. An N - N will still count as a mismatch
#[inline]
pub fn hamming_dist_none(alpha: &BaseBits, beta: &BaseBits) -> u32 {
    let nbits_and = alpha.nbits & beta.nbits;
    (((alpha.code ^ beta.code) & nbits_and).count_ones() / ENCODING_DIST)
        + ((!nbits_and).count_ones() / ENCODING_LENGTH)
}

/// Extract 'k' bits from the end of a u64 integer
#[inline]
fn extract_bits(n: u64, k: u32) -> u64 {
    !(!0u64 << k) & n
}

// Hamming distance functions that don't depend on BaseBits types
pub mod hamming {

    /// If you have encoded a string using hamming codes, this should work
    #[inline]
    pub fn hamming_code(alpha: u64, beta: u64) -> u32 {
        (alpha ^ beta).count_ones() / 2
    }

    /// Classical hamming distance on strings. Skips length check and will stop comparing after
    /// alpha is exhuasted
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
        assert_eq!(hamming_dist_nany(&alpha, &beta), 1);
    }

    #[test]
    fn test_single_bases() {
        // Selfs
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"A").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"G").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"C").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"T").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"N").unwrap(), &BaseBits::new(b"N").unwrap()),
            0
        );

        // Others
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"N").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"N").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"N").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_nany(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"N").unwrap()),
            0
        );
        // Selfs
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"A").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"G").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"C").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"T").unwrap()),
            0
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"N").unwrap(), &BaseBits::new(b"N").unwrap()),
            1
        );

        // Others
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"A").unwrap(), &BaseBits::new(b"N").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"G").unwrap(), &BaseBits::new(b"N").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"T").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"C").unwrap(), &BaseBits::new(b"N").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"A").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"G").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"C").unwrap()),
            1
        );
        assert_eq!(
            hamming_dist_none(&BaseBits::new(b"T").unwrap(), &BaseBits::new(b"N").unwrap()),
            1
        );
    }

    #[test]
    fn test_case_n_hamming() {
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"NCTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );
    }

    #[test]
    fn test_cases_bb_hamming() {
        // Test N encoding
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"ACTN").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );

        assert_eq!(
            hamming_dist_none(
                &BaseBits::new(b"ACTN").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            1
        );

        // Test N encoding
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"NCTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );

        assert_eq!(
            hamming_dist_none(
                &BaseBits::new(b"NCTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            1
        );

        // Test that unkown chars treated like Ns
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"ACT9").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );

        assert_eq!(
            hamming_dist_none(
                &BaseBits::new(b"ACT9").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            1
        );
        // Test regular equality
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"ACTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );
        assert_eq!(
            hamming_dist_none(
                &BaseBits::new(b"ACTG").unwrap(),
                &BaseBits::new(b"ACTG").unwrap()
            ),
            0
        );

        // Test Other string
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"GATACA").unwrap(),
                &BaseBits::new(b"GATACT").unwrap()
            ),
            1
        );
        assert_eq!(
            hamming_dist_none(
                &BaseBits::new(b"GATACA").unwrap(),
                &BaseBits::new(b"GATACT").unwrap()
            ),
            1
        );
        assert_eq!(
            hamming_dist_nany(
                &BaseBits::new(b"CATACAGATACTTCCATAGCT").unwrap(),
                &BaseBits::new(b"GATACAGATACTTCCATAGCA").unwrap()
            ),
            2
        );
        assert_eq!(
            hamming_dist_none(
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
        let beta = BaseBits::new(b"ACTGN").unwrap();
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
