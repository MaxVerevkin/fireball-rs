use super::common::CircleSize;
use crate::IntoScope;

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub(crate) size: CircleSize,
}

impl Sphere {
    pub fn new(size: impl Into<CircleSize>) -> Self {
        Self { size: size.into() }
    }
}

impl IntoScope for Sphere {
    fn into_scope(self, scope: &mut crate::Scope) {
        scope.0.push(crate::element::Element::Sphere(self));
    }
}
