use crate::Scope;
use std::fmt;
use std::ops::{Deref, DerefMut};

#[derive(Debug, Default)]
pub struct Document {
    root: Scope,
    detail: Option<u32>,
}

impl Document {
    pub fn set_detail(&mut self, detail: u32) {
        self.detail = Some(detail);
    }
}

impl Deref for Document {
    type Target = Scope;
    fn deref(&self) -> &Self::Target {
        &self.root
    }
}

impl DerefMut for Document {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.root
    }
}

impl fmt::Display for Document {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(detail) = self.detail {
            writeln!(f, "$fn = {detail};")?;
        }
        for obj in &self.root.0 {
            writeln!(f, "{obj}")?;
        }
        Ok(())
    }
}
