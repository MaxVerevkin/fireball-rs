## Brief
This program aims to reconstruct the trajectory of a fireball given
only data collected from witnesses.

Run:
```
$ fireball -h
```
to see which parameters of the implementation can be tweaked.

## Implementation
The current implementation uses binary search to find two points
in the space that represents the trajectory of a fireball. The binary search
itself minimizes the mean-square-error of observations given trajectory.

## Data format
The data is given in a form:
```
    lat lon h A z_start h_start z_end h_end t
```
where:
* `lat` - geographical latitude of observer's location,
* `lon` - geographical longitude of observer's location,
* `h` - the observer's altitude above the sea,
* `A` - observed descent angle,
* `z_start` - azimuth of the begining of the observation,
* `h_start` - altitude of the begining of the observation,
* `z_end` - azimuth of the end of the observation,
* `h_end` - altitude of the end of the observation,
* `t` - the duration of the observation.

The file with data should contain only the numbers, one observation per line. The example of such a file is prvided as `data-2020-01-21` (data from a real event)
