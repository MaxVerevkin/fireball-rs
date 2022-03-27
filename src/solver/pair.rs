use crate::data::DataSample;
use crate::structs::Line;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairTrajectory {
    pub line: Line,
    pub weight: f64,
}

impl PairTrajectory {
    pub fn calculate(s1: &DataSample, s2: &DataSample) -> Option<Self> {
        // Skip this pair if their obsercations are too close
        //
        // Note: Ten meters is a randomly chosen number, there is no deep meaning to it. Perhaps
        // this check should be replaced with something more sophisticated.
        if (s1.location - s2.location).len() < 10.0 {
            return None;
        }

        let axis1 = s1.axis?;
        let plane1 = s1.plane?;
        let plane2 = s2.plane?;

        let mut dir = plane1.cross(plane2);

        // weight = sin^2(plane1^plane2)
        //
        // Note: Perpendicular planes are expected to give more reliable results, while almost
        // parallel planes should be ignored.
        let weight = dir.len_sq();

        dir.normalize();
        if !dir.is_normal() {
            // Planes are parallel
            return None;
        }

        // Check direction
        let dir = match (s1.direction_matches(dir), s2.direction_matches(dir)) {
            // Observers are giving opposite directions
            (true, false) | (false, true) => return None,
            // Both observers argee that velocity should be negated
            (false, false) => -dir,
            // Both observers argee that velocity should not be changed
            (true, true) => dir,
        };

        let l1 = (s2.location - s1.location).dot(plane2) / axis1.dot(plane2);
        let point = s1.location + axis1 * l1;

        // Ensure that the point actually belongs to both planes
        assert_approx_eq!(plane1.dot((point - s1.location).normalized()), 0.0);
        assert_approx_eq!(plane2.dot((point - s2.location).normalized()), 0.0);

        (s1.observation_matches(point, dir) && s2.observation_matches(point, dir)).then(|| Self {
            line: Line {
                point,
                direction: dir,
            },
            weight,
        })
    }
}
