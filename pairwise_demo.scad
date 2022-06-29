///

p1 = [10, 40, 0];
begin1_az = [130, 60];
end1_az = [180, 40];

p2 = [20, -30, 0];
begin2_az = [10, 60];
end2_az = [-30, 40];

///

translate(p1) sphere(2, $fn = 55);
translate(p2) sphere(2, $fn = 55);

begin1 = az_to_xyz(begin1_az);
end1 = az_to_xyz(end1_az);
begin2 = az_to_xyz(begin2_az);
end2 = az_to_xyz(end2_az);

plane1 = normalize(cross(end1, begin1));
plane2 = normalize(cross(az_to_xyz(begin2_az), az_to_xyz(end2_az)));
dir = normalize(cross(plane1, plane2));

l_begin1 = dot(p2 - p1, plane2) / dot(begin1, plane2);
l_end1 = dot(p2 - p1, plane2) / dot(end1, plane2);
l_begin2 = dot(p1 - p2, plane1) / dot(begin2, plane1);
l_end2 = dot(p1 - p2, plane1) / dot(end2, plane1);

point = (p1 + end1 * l_end1 + p2 + end2 * l_end2) / 2;

arrow(p1, begin1_az, l_begin1);
arrow(p1, end1_az, l_end1);
arrow(p2, begin2_az, l_begin2);
arrow(p2, end2_az, l_end2);

color("green") {
    translate(point) sphere(1, $fn=55);

    hull() {
        translate(point) sphere(0.4, $fn=55);
        translate(point + dir * 1000) sphere(0.4, $fn=55);
    }
}

module arrow(from, direction, size) {
    translate(from)
    rotate([0, 0, -direction[0]])
    rotate([direction[1]-90, 0, 0])
    {
        translate([0, 0, size - 3]) cylinder(r1 = 2, r2 = 0, h = 3, $fn = 20);
        translate([-1 / 2, -1 / 2, 0]) cube([1, 1, size - 3]);
    }
}

function az_to_xyz(az) = [sin(az[0])*cos(az[1]), cos(az[0])*cos(az[1]), sin(az[1])];

function normalize(v) = v / sqrt(dot(v, v));

function dot(u, v) = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
