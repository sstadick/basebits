#[cfg(test)]
mod tests {
    use super::hamming::*;
    use super::encode_dna::*;
    #[test]
    fn test_hamming_dist() {
        assert_eq!(hamming_str("ACTG", "ACTT"), 1);
        assert_eq!(hamming_str("ACTG", "ACTTT"), 1);
    }

    #[test]
    fn test_hamming_code() {
        assert_eq!(hamming_code(encode("ACTG"), encode("ACTT")), 1);
    }

    #[test]
    fn test_encode() {
        let code = encode("ACTG");
        println!("{:b}", code);
        assert_eq!(code, 0b000110101011);
    }
}

// Encode

pub mod encode_dna {
    pub fn encode(seq: &str) -> u64 {
        let mut code: u64 = 0;
        for c in seq.chars() {
            code = (code << 3) | match c {
                'A' => 0b000,
                'C' => 0b110,
                'T' => 0b101,
                'G' => 0b011,
                'N' => 0b100,
                 _ => 0b111,
            }
        }
        code
    }
}

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
