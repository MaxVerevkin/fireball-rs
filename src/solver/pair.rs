use crate::approx_eq::ApproxEq;

use common::nalgebra::Unit;
use common::obs_data::DataSample;
use common::structs::Line;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairTrajectory {
    /// The calculated trajectory
    pub line: Line,
    /// How trustworthy this trajectory is
    pub weight: f64,
}

impl PairTrajectory {
    pub fn calculate(s1: &DataSample, s2: &DataSample) -> Option<Self> {
        // Skip this pair if their obsercations are too close
        //
        // Note: Ten meters is a randomly chosen number, there is no deep meaning to it. Perhaps
        // this check should be replaced with something more sophisticated.
        if (s1.location - s2.location).norm() < 10.0 {
            return None;
        }

        let axis1 = s1.axis?;
        let plane1 = s1.plane?;
        let plane2 = s2.plane?;

        // weight = sin^2(plane1^plane2)
        //
        // Note: Perpendicular planes are expected to give more reliable results, while almost
        // parallel planes should be ignored.
        let (direction, mut weight) = Unit::new_and_get(plane1.cross(&plane2));
        weight *= weight;

        if weight.approx_eq(0.0) {
            // Planes are parallel
            return None;
        }

        // Check direction
        let direction = match (
            s1.direction_matches(direction),
            s2.direction_matches(direction),
        ) {
            // Observers are giving opposite directions
            (true, false) | (false, true) => return None,
            // Both observers argee that velocity should be negated
            (false, false) => -direction,
            // Both observers argee that velocity should not be changed
            (true, true) => direction,
        };

        let l1 = (s2.location - s1.location).dot(&plane2) / axis1.dot(&plane2);
        let point = s1.location + axis1.into_inner() * l1;

        // Ensure that the point actually belongs to both planes
        assert_approx_eq!(plane1.dot(&(point - s1.location).normalize()), 0.0);
        assert_approx_eq!(plane2.dot(&(point - s2.location).normalize()), 0.0);

        if !s1.observation_matches(Line { point, direction })
            || !s2.observation_matches(Line { point, direction })
        {
            return None;
        }

        Some(Self {
            line: Line { point, direction },
            weight,
        })
    }
}
