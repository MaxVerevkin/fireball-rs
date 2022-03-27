use crate::element::Element;
use crate::objects::*;
use crate::Vec3;
use std::fmt;

#[derive(Debug, Default)]
pub struct Scope(pub(crate) Vec<Element>);

macro_rules! get_last {
    ($self:ident, $t:ident) => {
        match $self.0.last_mut().unwrap() {
            Element::$t(var) => var.borrow(),
            _ => unreachable!(),
        }
    };
}

macro_rules! new_scope {
    ($into:ident) => {{
        let mut scope = Scope::default();
        $into.into_scope(&mut scope);
        scope
    }};
}

impl Scope {
    pub fn sphere(&mut self, size: impl Into<CircleSize>) {
        self.0.push(Element::Sphere(Sphere::new(size)));
    }

    pub fn cube(&mut self, size: impl Into<CubeSize>) -> Cube<'_> {
        self.0.push(Element::Cube(Cube::new(size)));
        get_last!(self, Cube)
    }

    pub fn cylinder(&mut self, height: f32, size: impl Into<CircleSize>) -> Cylinder<'_> {
        self.0.push(Element::Cylinder(Cylinder::new(height, size)));
        get_last!(self, Cylinder)
    }

    pub fn translate(&mut self, arg: impl Into<Vec3>, scope: impl IntoScope) {
        self.0
            .push(Element::Translate(arg.into(), new_scope!(scope)));
    }

    pub fn rotate(&mut self, arg: impl Into<Vec3>, scope: impl FnOnce(&mut Scope)) {
        self.0.push(Element::Rotate(arg.into(), new_scope!(scope)));
    }

    pub fn color(&mut self, color: impl Into<Color>, alpha: Option<f32>, scope: impl IntoScope) {
        self.0
            .push(Element::Color(color.into(), alpha, new_scope!(scope)));
    }

    pub fn hull(&mut self, scope_fn: impl FnOnce(&mut Scope)) {
        let mut scope = Scope::default();
        scope_fn(&mut scope);
        self.0.push(Element::Hull(scope));
    }
}

impl fmt::Display for Scope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0.as_slice() {
            [] => f.write_str("{}\n"),
            [obj] => writeln!(f, "{obj}"),
            objects => {
                f.write_str("{\n")?;
                for obj in objects {
                    writeln!(f, "{obj}")?;
                }
                f.write_str("}\n")
            }
        }
    }
}

pub trait IntoScope {
    fn into_scope(self, scope: &mut Scope);
}

impl<T: FnOnce(&mut Scope)> IntoScope for T {
    fn into_scope(self, scope: &mut Scope) {
        self(scope)
    }
}

#[macro_export]
macro_rules! scope {
    ($ctx:ident, $code:expr) => {
        |$ctx: &mut $crate::Scope| $code
    };
}
