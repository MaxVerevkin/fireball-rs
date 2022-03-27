use crate::scad_value::ScadValue;

#[derive(Debug, Clone, Copy)]
pub enum CircleSize {
    Radius(f32),
    Diameter(f32),
}

impl From<f32> for CircleSize {
    fn from(r: f32) -> Self {
        Self::Radius(r)
    }
}

impl ScadValue for CircleSize {
    fn scad_str(&self) -> String {
        match *self {
            Self::Radius(r) => format!("r={r}"),
            Self::Diameter(d) => format!("d={d}"),
        }
    }
}
