extern crate hammy;
use hammy::{encode, dist};

fn main() {
    let string1 = "ACTGACTG".to_string();
    let string2 = "ACTTACTG".to_string();

    let string1 = encode(string1);
    let string2 = encode(string2);

    assert_eq!(dist(string1, string2), 1);
}
