//! Useful constants to have

/// The radius of the Earth in meters
///
/// More prescription is useless
/// since the radius can vary.
pub const EARTH_R: f64 = 6_371_000.;

/// Unicode's degree symbol
pub const DEGREE_SYM: char = '\u{00b0}';

/// Unicode's delta symbol
pub const DELTA_SYM: char = '\u{0394}';

/// Unicode's sigma symbol
pub const SIGMA_SYM: char = '\u{03c3}';

/// Unicode's plus-minus symbol
pub const PLUS_MINUS: char = '\u{00b1}';

/// Go to the begining of the previous line
pub const ANSI_GOTO_PREV_LINE: &str = "\u{001b}[1F";

/// Clear current line
pub const ANSI_CLEAR_LINE: &str = "\u{001b}[K";
