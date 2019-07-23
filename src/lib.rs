#[cfg(test)]
mod tests {
    use super::hamming::*;
    use super::base_bits::*;
    #[test]
    fn test_hamming_dist() {
        assert_eq!(hamming_str("ACTG", "ACTT"), 1);
        assert_eq!(hamming_str("ACTG", "ACTTT"), 1);
    }

    #[test]
    fn test_base_bits() {
        let alpha = BaseBits::new("ACTG");
        let beta = BaseBits::new("ACTT");
        assert_eq!(hamming_dist(&alpha, &beta), 1);
    }

    #[test]
    fn test_encoding() {
        let bb = BaseBits::new("ACTG");
        assert_eq!(bb.code, 0b000110101011);
    }
}

/// Encode a DNA string of up to 21 bases as a u64 for fast hamming distance calculations.
/// TODO: Add a bump to use u128 or maybe bigint if 21 chars is not enough. 
pub mod base_bits {
    pub const ENCODING_DIST: u32 = 2;
    pub const ENCODING_LENGTH: u32 = 3;
    pub const UNDETERMINED: u64 = 0b100;
    pub const ANY: u64 = 0b111;

    struct Bases;
    impl Bases {
        const A: u64 = 0b000;
        const C: u64 = 0b110;
        const T: u64 = 0b101;
        const G: u64 = 0b011;
        const N: u64 = UNDETERMINED;
        const STAR: u64 = ANY;
    }
  
    //#[derive(Copy, Clone)]
    pub struct BaseBits {
        pub code: u64
    }

    impl BaseBits {
        pub fn new(seq: &str)-> BaseBits {
            let mut code: u64 = 0;
            for c in seq.chars() {
                code = (code << ENCODING_LENGTH) | match c {
                    'A' => Bases::A,
                    'C' => Bases::C,
                    'T' => Bases::T,
                    'G' => Bases::G,
                    'N' => Bases::N,
                     _ => Bases::STAR,
                }
            }
            BaseBits{code}
        }
    }

    #[inline]
    pub fn hamming_dist(alpha: &BaseBits, beta: &BaseBits) -> u32 {
        (alpha.code ^ beta.code).count_ones() / 2
    }
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
