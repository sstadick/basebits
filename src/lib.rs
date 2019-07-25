#[cfg(test)]
mod tests {
    use super::hamming::*;
    use super::base_bits::*;
    #[test]
    fn test_hamming_str_dist() {
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
    fn test_cases_bb_hamming() {
        // Test N encoding
        assert_eq!(hamming_dist(&BaseBits::new("ACTN"), &BaseBits::new("ACTG")), 1);
        // Test * encoding
        assert_eq!(hamming_dist(&BaseBits::new("ACT*"), &BaseBits::new("ACTG")), 0);
        // Test that unkown chars treated like Ns
        assert_eq!(hamming_dist(&BaseBits::new("ACT9"), &BaseBits::new("ACTG")), 1);
        // Test regular equality
        assert_eq!(hamming_dist(&BaseBits::new("ACTG"), &BaseBits::new("ACTG")), 0);

        // Test Other string
        assert_eq!(hamming_dist(&BaseBits::new("GATACA"), &BaseBits::new("GATACT")), 1);
    }

    #[test]
    fn test_bb_to_string() {
        let alpha = BaseBits::new("GCTAN");
        let beta = BaseBits::new("ACTG*");
        println!("Alpha: {}", alpha);
        println!("Beta: {}", beta);
        assert_eq!(alpha.to_string(), "GCTAN".to_string());
        assert_eq!(beta.to_string(), "ACTG*".to_string());
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
    use std::fmt;

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
        pub code: u64,
        len: usize
    }

    impl BaseBits {
        pub fn new(seq: &str)-> BaseBits {
            let mut code: u64 = 0;
            let len = seq.len();
            for c in seq.chars() {
                code = (code << ENCODING_LENGTH) | match c {
                    'A' => Bases::A,
                    'C' => Bases::C,
                    'T' => Bases::T,
                    'G' => Bases::G,
                    'N' => Bases::N,
                    '*' => Bases::STAR,
                     _ => Bases::N,
                }
            }
            BaseBits{code, len}
        }

        fn decode(&self) -> String {
            // firgure out how to pop off ENCODING_Length bits at a time
            // and decode those bits... 
            let mut s = String::from("");
            let mut code = self.code;
            println!("Code: {:#b}", code);
            for _ in 0..self.len {
                let base = extract_bits(code, ENCODING_LENGTH);
                println!("Base: {:#b}", base);
                code = code >> ENCODING_LENGTH;
                println!("Code: {:#b}", code);
                s.push(match base {
                    Bases::A => 'A',
                    Bases::C => 'C',
                    Bases::T => 'T',
                    Bases::G => 'G',
                    Bases::N => 'N',
                    Bases::STAR => '*',
                    _ => 'N'
                });
            }
            s.chars().rev().collect()
        }
    }

    impl fmt::Display for BaseBits {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}", self.decode())
        }
    }

    #[inline]
    pub fn hamming_dist(alpha: &BaseBits, beta: &BaseBits) -> u32 {
        (alpha.code ^ beta.code).count_ones() / 2
    }

    // Util functions
    /// Extract 'k' bits from the end of a u64 integer
    #[inline]
    fn extract_bits(n: u64, k: u32) -> u64 {
        let mut extractor: u64 = 0;
        for _ in 0..k {
            extractor <<= 1;
            extractor |= 1;
        }
        extractor & n 
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
