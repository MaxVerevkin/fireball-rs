use crate::Vec3;

pub trait ScadValue {
    fn scad_str(&self) -> String;
}

impl ScadValue for Vec3 {
    fn scad_str(&self) -> String {
        format!("[{},{},{}]", self.x, self.y, self.z)
    }
}
