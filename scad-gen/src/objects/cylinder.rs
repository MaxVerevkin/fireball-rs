use super::common::CircleSize;

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct CylinderInner {
    pub(crate) height: f32,
    pub(crate) size: CircleSize,
    pub(crate) center: bool,
}

gen_wrapper_type!(Cylinder, CylinderInner);

impl Cylinder<'static> {
    pub fn new(height: f32, size: impl Into<CircleSize>) -> Self {
        Self::Owned(CylinderInner {
            height,
            size: size.into(),
            center: false,
        })
    }
}

impl Cylinder<'_> {
    pub fn center(mut self) -> Self {
        self.center = true;
        self
    }
}
