use crate::objects::*;
use crate::scad_value::ScadValue;
use crate::Scope;
use crate::Vec3;
use std::fmt;

#[derive(Debug)]
pub enum Element {
    // 3D
    Sphere(Sphere),
    Cube(Cube<'static>),
    Cylinder(Cylinder<'static>),
    // Transformations
    Translate(Vec3, Scope),
    Rotate(Vec3, Scope),
    Color(Color, Option<f32>, Scope),
    Hull(Scope),
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Sphere(sphere) => write!(f, "sphere({});", sphere.size.scad_str()),
            Self::Cube(cube) => {
                write!(f, "cube({},center={});", cube.size.scad_str(), cube.center)
            }
            Self::Cylinder(cylinder) => {
                write!(
                    f,
                    "cylinder(h={},{},center={});",
                    cylinder.height,
                    cylinder.size.scad_str(),
                    cylinder.center
                )
            }
            Self::Translate(arg, scope) => write!(f, "translate({}) {scope}", arg.scad_str()),
            Self::Rotate(arg, scope) => write!(f, "rotate({}) {scope}", arg.scad_str()),
            Self::Color(color, None, scope) => {
                write!(f, "color({}) {scope}", color.scad_str())
            }
            Self::Color(color, Some(alpha), scope) => {
                write!(f, "color({},{alpha}) {scope}", color.scad_str())
            }
            Self::Hull(scope) => write!(f, "hull() {scope}"),
        }
    }
}
