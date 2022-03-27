use crate::objects::*;
use crate::Vec3;

#[inline(always)]
pub fn vec3(x: f32, y: f32, z: f32) -> Vec3 {
    Vec3::new(x, y, z)
}

#[inline(always)]
pub fn cube(size: impl Into<CubeSize>) -> Cube<'static> {
    Cube::new(size)
}

#[inline(always)]
pub fn sphere(size: impl Into<CircleSize>) -> Sphere {
    Sphere::new(size)
}

#[inline(always)]
pub fn cylinder(height: f32, size: impl Into<CircleSize>) -> Cylinder<'static> {
    Cylinder::new(height, size)
}
