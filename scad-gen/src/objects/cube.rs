use crate::scad_value::ScadValue;
use crate::Vec3;

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct CubeInner {
    pub(crate) size: CubeSize,
    pub(crate) center: bool,
}

gen_wrapper_type!(Cube, CubeInner);

impl Cube<'static> {
    pub fn new(size: impl Into<CubeSize>) -> Self {
        Self::Owned(CubeInner {
            size: size.into(),
            center: false,
        })
    }
}

impl Cube<'_> {
    pub fn center(mut self) -> Self {
        self.center = true;
        self
    }
}

#[derive(Debug, Clone, Copy)]
pub enum CubeSize {
    Size(f32),
    Dimensions(Vec3),
}

impl From<f32> for CubeSize {
    fn from(size: f32) -> Self {
        Self::Size(size)
    }
}

impl<T: Into<Vec3>> From<T> for CubeSize {
    fn from(dims: T) -> Self {
        Self::Dimensions(dims.into())
    }
}

impl ScadValue for CubeSize {
    fn scad_str(&self) -> String {
        match *self {
            Self::Size(size) => size.to_string(),
            Self::Dimensions(dims) => dims.scad_str(),
        }
    }
}
