use common::obs_data::{DataSample, PartDiff};
use common::structs::{Geodetic, Line, UnitVec3, Vec3};

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
        let axis2 = s2.axis?;
        let plane1 = s1.plane?;
        let plane2 = s2.plane?;

        let dot1 = axis1.dot(*plane2).abs();
        let dot2 = axis2.dot(*plane1).abs();

        const MIN_ACCEPTABLE_SIN: f64 = 0.1; // ~ sin(6 deg)
        if dot1 < MIN_ACCEPTABLE_SIN || dot2 < MIN_ACCEPTABLE_SIN {
            return None;
        }

        let weight = (dot1 * dot2).sqrt();

        // Check direction
        let direction = UnitVec3::new_normalize(plane1.cross(*plane2));
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

        let l1 = (s2.location - s1.location).dot(*plane2) / axis1.dot(*plane2);
        let point1 = s1.location + axis1 * l1;
        let l2 = (s1.location - s2.location).dot(*plane1) / axis2.dot(*plane1);
        let point2 = s2.location + axis2 * l2;
        let point = (point1 + point2) * 0.5;

        if axis1.dot(*axis2).acos() < 10f64.to_radians() {
            if l1 > 5e6 && l2 > 5e6 {
                dbg!(
                    weight,
                    l1,
                    l2,
                    Geodetic::from_geocentric_cartesian(point, 10)
                );
            }
        }
        // Ensure that the point actually belongs to both planes
        assert_approx_eq!(plane1.dot((point - s1.location).normalize()), 0.0);
        assert_approx_eq!(plane2.dot((point - s2.location).normalize()), 0.0);

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

    /// Calculate the derivative of the `calculate` function is respect to `s1`'s `pd`.
    ///
    /// Call this function only if `calculate` have returned `Some`.
    ///
    /// Returns (point_diff, direction_diff, weight_diff)
    pub fn diff(s1: &DataSample, s2: &DataSample, pd: &PartDiff) -> Option<(Vec3, Vec3, f64)> {
        let axis1 = s1.axis?;
        let axis2 = s2.axis?;
        let plane1 = s1.plane?;
        let plane2 = s2.plane?;

        let axis1_diff = s1.axis_diff(pd)?;
        let plane1_diff = s1.plane_diff(pd)?;

        let dot1 = axis1.dot(*plane2).abs();
        let dot1_diff = axis1_diff.dot(*plane2).abs();
        let dot2 = axis2.dot(*plane1).abs();
        let dot2_diff = axis2.dot(plane1_diff).abs();

        let weight_sq = dot1 * dot2;
        let weight_sq_diff = dot1 * dot2_diff + dot1_diff * dot2;
        let weight = weight_sq.sqrt();
        let weight_diff = 0.5 / weight * weight_sq_diff;

        let cross = plane1.cross(*plane2);
        let cross_diff = plane1_diff.cross(*plane2);

        let direction = UnitVec3::new_normalize(cross);
        let direction_diff = cross.normalize_diff(cross_diff);

        // Check direction
        let (_direction, direction_diff) = match (
            s1.direction_matches(direction),
            s2.direction_matches(direction),
        ) {
            // Observers are giving opposite directions
            (true, false) | (false, true) => return None,
            // Both observers argee that velocity should be negated
            (false, false) => (-direction, -direction_diff),
            // Both observers argee that velocity should not be changed
            (true, true) => (direction, direction_diff),
        };

        let l1 = (s2.location - s1.location).dot(*plane2) / axis1.dot(*plane2);
        let l1_diff = -(s2.location - s1.location).dot(*plane2) / axis1.dot(*plane2).powi(2)
            * axis1_diff.dot(*plane2);

        let l2 = (s1.location - s2.location).dot(*plane1) / axis2.dot(*plane1);
        let l2_diff = (s1.location - s2.location).dot(plane1_diff) / axis2.dot(*plane1)
            - (s1.location - s2.location).dot(*plane1) / axis2.dot(*plane1).powi(2)
                * axis2.dot(plane1_diff);

        let point1 = s1.location + axis1 * l1;
        let point1_diff = axis1 * l1_diff + axis1_diff * l1;

        let point2 = s2.location + axis2 * l2;
        let point2_diff = axis2 * l2_diff;

        let _point = (point1 + point2) * 0.5;
        let point_diff = (point1_diff + point2_diff) * 0.5;

        Some((point_diff, direction_diff, weight_diff))
    }
}
