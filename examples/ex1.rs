extern crate basebits;
use basebits::{hamming_dist_nany, BaseBits};

fn main() {
    let string1 = b"ACTGACTG";
    let string2 = b"ACTTACTG";

    let string1 = BaseBits::new(string1).unwrap();
    let string2 = BaseBits::new(string2).unwrap();

    assert_eq!(hamming_dist_nany(&string1, &string2), 1);
}
