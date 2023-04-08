use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fireball::maths::{lambda, lambda_old};
use fireball::structs::Vec3;

fn old(c: &mut Criterion) {
    c.bench_function("old", |b| {
        b.iter(|| {
            lambda_old(
                black_box(Vec3::default()),
                black_box(Vec3::x_axis()),
                black_box(Vec3::x_axis().into_inner() + Vec3::y_axis().into_inner()),
                black_box(Vec3::z_axis()),
            )
        })
    });
}

fn new(c: &mut Criterion) {
    c.bench_function("new", |b| {
        b.iter(|| {
            lambda(
                black_box(Vec3::default()),
                black_box(Vec3::x_axis()),
                black_box(Vec3::x_axis().into_inner() + Vec3::y_axis().into_inner()),
                black_box(Vec3::z_axis()),
            )
        })
    });
}

criterion_group!(benches, old, new);
criterion_main!(benches);
