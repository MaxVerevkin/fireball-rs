pub mod constructors;
mod document;
mod element;
pub mod objects;
mod scad_value;
mod scope;
mod vector;

pub use document::*;
pub use objects::CircleSize::*;
pub use scope::*;
pub use vector::*;

#[cfg(test)]
mod tests {
    use super::constructors::*;
    use super::*;

    #[test]
    fn simple_cube() {
        let mut doc = Document::default();
        doc.cube(1f32).center();
        doc.cube(vec3(0.0, 1.0, 3.0));
        assert_eq!(
            doc.to_string(),
            "cube(1,center=true);\ncube([0,1,3],center=false);\n"
        );
    }

    #[test]
    fn into_scope() {
        let mut doc1 = Document::default();
        doc1.translate(vec3(1., 0., -1.), |ctx: &mut Scope| {
            ctx.cube(1f32).center();
        });

        let mut doc2 = Document::default();
        doc2.translate(vec3(1., 0., -1.), cube(1f32).center());

        assert_eq!(doc1.to_string(), doc2.to_string());
    }

    #[test]
    fn scope_macro() {
        let mut doc1 = Document::default();
        doc1.color("green", None, |ctx: &mut Scope| {
            ctx.cube(1f32);
            ctx.translate(vec3(1., 1., 1.), cube(1f32));
        });

        let mut doc2 = Document::default();
        doc2.color(
            "green",
            None,
            scope!(ctx, {
                ctx.cube(1f32);
                ctx.translate(vec3(1., 1., 1.), cube(1f32));
            }),
        );

        assert_eq!(doc1.to_string(), doc2.to_string());
    }
}
