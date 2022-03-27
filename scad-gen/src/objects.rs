use crate::scad_value::ScadValue;
use crate::Vec3;
use std::borrow::Cow;

macro_rules! gen_wrapper_type {
    ($outer_type:ident, $inner_type:ident) => {
        #[derive(Debug)]
        pub enum $outer_type<'a> {
            Owned($inner_type),
            Borrowed(&'a mut $inner_type),
        }

        impl ::std::ops::Deref for $outer_type<'_> {
            type Target = $inner_type;
            fn deref(&self) -> &Self::Target {
                match self {
                    Self::Owned(inner) => inner,
                    Self::Borrowed(borrowed) => borrowed,
                }
            }
        }

        impl ::std::ops::DerefMut for $outer_type<'_> {
            fn deref_mut(&mut self) -> &mut Self::Target {
                match self {
                    Self::Owned(inner) => inner,
                    Self::Borrowed(borrowed) => borrowed,
                }
            }
        }

        impl $outer_type<'static> {
            pub(crate) fn borrow(&mut self) -> $outer_type<'_> {
                $outer_type::Borrowed(self)
            }
        }

        impl $crate::IntoScope for $outer_type<'static> {
            fn into_scope(self, scope: &mut $crate::Scope) {
                scope.0.push($crate::element::Element::$outer_type(self));
            }
        }
    };
}

mod common;
mod cube;
mod cylinder;
mod sphere;

pub use common::*;
pub use cube::*;
pub use cylinder::*;
pub use sphere::*;

#[derive(Debug, Clone)]
pub enum Color {
    Rgb(Vec3),
    String(Cow<'static, str>),
}

impl<T: Into<Vec3>> From<T> for Color {
    fn from(rgb: T) -> Self {
        Self::Rgb(rgb.into())
    }
}

impl From<&'static str> for Color {
    fn from(text: &'static str) -> Self {
        Self::String(Cow::Borrowed(text))
    }
}

impl From<String> for Color {
    fn from(text: String) -> Self {
        Self::String(Cow::Owned(text))
    }
}

impl ScadValue for Color {
    fn scad_str(&self) -> String {
        match self {
            Color::Rgb(rgb) => rgb.scad_str(),
            Color::String(text) => format!("\"{text}\""),
        }
    }
}
