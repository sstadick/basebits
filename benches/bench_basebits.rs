#[macro_use]
extern crate criterion; 

use criterion::Criterion;
use criterion::black_box;

extern crate basebits;
use basebits::*;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("ham distance classical on strings",
                     |b| b.iter(||
                      hamming::hamming_str(black_box("ACTGACTGACTG"),
                                           black_box("ACTGGGGGACTG"))));
    let a = BaseBits::new(b"ACTGACTGACTG").unwrap();
    let be = BaseBits::new(b"ACTGGGGGACTG").unwrap();
    c.bench_function("BaseBits pre encoded input",
                      move |b| b.iter(||
                      hamming_dist(&a, &be)));
    
    c.bench_function("BaseBits with encoding",
                     |b| b.iter(||
                      hamming_dist(&BaseBits::new(b"ACTGACTGACTG").unwrap(), &BaseBits::new(b"ACTGGGGGACTG").unwrap())));
    c.bench_function("BaseBits encoding cost",
                     |b| b.iter(||
                      black_box(BaseBits::new(b"ACTGACTGACTG").unwrap())));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
