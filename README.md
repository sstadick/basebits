![docs](https://docs.rs/basebits/badge.svg)

# basebits

A library for memory efficient short DNA sequence encoding.

## Synopsis

When to use this library? 
If you are comparing strings against each other more than 4 times, it
becomes more efficient to pay the cost of encoding them. 

## Operations

Constant time hamming distance calculations.

## Example

```rust
use basebits::{BaseBits, hamming_dist};

fn main() {
    let string1 = b"ACTGACTG";
    let string2 = b"ACTTACTG";

    let string1 = BaseBits::new(string1).unwrap();
    let string2 = BaseBits::new(string2).unwrap();

    assert_eq!(hamming_dist(&string1, &string2), 1);
}
```

## Reference

See 'Constant Time Hamming Distance' section:
https://www.biorxiv.org/content/10.1101/648683v1.full

## Future directions

FFT stuff?

