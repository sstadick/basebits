#[macro_use]
extern crate criterion; 

use criterion::Criterion;
use criterion::black_box;

extern crate basebits;
use basebits::*;

fn criterion_benchmark(c: &mut Criterion) {
    let alpha = BaseBits::new("ACTGACTGACTG").unwrap();
    let beta = BaseBits::new("ACTGGGGGACTG").unwrap();
    c.bench_function("ham xor pre encoded",
                     move |b| b.iter(||
                      hamming::hamming_code(alpha.code, beta.code)));
    c.bench_function("ham xor",
                     |b| b.iter(||
                      hamming::hamming_code(BaseBits::new("ACTGACTGACTG").unwrap().code,
                      BaseBits::new("ACTGGGGGACTG").unwrap().code)));
    c.bench_function("ham str",
                     |b| b.iter(||
                      hamming::hamming_str(black_box("ACTGACTGACTG"),
                                           black_box("ACTGGGGGACTG"))));
    let a = BaseBits::new("ACTGACTGACTG").unwrap();
    let be = BaseBits::new("ACTGGGGGACTG").unwrap();
    c.bench_function("BaseBits Pre",
                      move |b| b.iter(||
                      hamming_dist(&a, &be)));
    
    c.bench_function("BaseBits with encoding",
                     |b| b.iter(||
                      hamming_dist(&BaseBits::new("ACTGACTGACTG").unwrap(), &BaseBits::new("ACTGGGGGACTG").unwrap())));
    c.bench_function("Encoding Cost",
                     |b| b.iter(||
                      black_box(BaseBits::new("ACTGACTGACTG").unwrap())));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
