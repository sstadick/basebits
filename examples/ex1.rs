extern crate basebits;
use basebits::{BaseBits, hamming_dist};

fn main() {
    let string1 = "ACTGACTG";
    let string2 = "ACTTACTG";

    let string1 = BaseBits::new(string1).unwrap();
    let string2 = BaseBits::new(string2).unwrap();

    assert_eq!(hamming_dist(&string1, &string2), 1);
}
