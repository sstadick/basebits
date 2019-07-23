#[macro_use]
extern crate criterion; 

use criterion::Criterion;
use criterion::black_box;

extern crate basebits;
use basebits::{hamming, base_bits::*};

fn criterion_benchmark(c: &mut Criterion) {
    let alpha = BaseBits::new("ACTGACTGACTG");
    let beta = BaseBits::new("ACTGGGGGACTG");
    c.bench_function("ham xor pre encoded",
                     move |b| b.iter(||
                      hamming::hamming_code(alpha.code, beta.code)));
    c.bench_function("ham xor",
                     |b| b.iter(||
                      hamming::hamming_code(BaseBits::new("ACTGACTGACTG").code,
                      BaseBits::new("ACTGGGGGACTG").code)));
    c.bench_function("ham str",
                     |b| b.iter(||
                      hamming::hamming_str(black_box("ACTGACTGACTG"),
                                           black_box("ACTGGGGGACTG"))));
    let a = BaseBits::new("ACTGACTGACTG");
    let be = BaseBits::new("ACTGGGGGACTG");
    c.bench_function("BaseBits Pre",
                      move |b| b.iter(||
                      hamming_dist(&a, &be)));
    
    c.bench_function("BaseBits with encoding",
                     |b| b.iter(||
                      hamming_dist(&BaseBits::new("ACTGACTGACTG"), &BaseBits::new("ACTGGGGGACTG"))));
    c.bench_function("Encoding Cost",
                     |b| b.iter(||
                      black_box(BaseBits::new("ACTGACTGACTG"))));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
