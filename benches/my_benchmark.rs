#[macro_use]
extern crate criterion; 

use criterion::Criterion;
use criterion::black_box;

extern crate hammy;
use hammy::{hamming, encode_dna};

fn criterion_benchmark(c: &mut Criterion) {
    let alpha = encode_dna::encode("ACTGACTGACTG");
    let beta = encode_dna::encode("ACTGGGGGACTG");
    c.bench_function("ham xor pre encoded",
                     move |b| b.iter(||
                      hamming::hamming_code(alpha, beta)));
    c.bench_function("ham xor",
                     |b| b.iter(||
                      hamming::hamming_code(encode_dna::encode("ACTGACTGACTG"),
                      encode_dna::encode("ACTGGGGGACTG"))));
    c.bench_function("ham str",
                     |b| b.iter(||
                      hamming::hamming_str(black_box("ACTGACTGACTG"),
                                           black_box("ACTGGGGGACTG"))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
