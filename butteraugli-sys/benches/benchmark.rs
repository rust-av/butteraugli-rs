use butteraugli_sys::butteraugli;
use criterion::{criterion_group, criterion_main, Criterion};

use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("butteraugli", |b| {
        b.iter(|| {
            butteraugli(
                "test-images/1a.png".to_string(),
                "test-images/1b.png".to_string(),
            )
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .measurement_time(Duration::new(8, 0))
        .sample_size(10);
    targets = criterion_benchmark
}

criterion_main!(benches);
