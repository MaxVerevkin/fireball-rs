## Brief
This program aims to reconstruct the trajectory of a fireball given
only data collected from witnesses.

Run:
```
$ fireball -h
```
to see which parameters of the implementation can be tweaked.

## Project Structure

- /src/maths.rs & /src/aprox_eq.rs & /src/constants.rs - Some mathematical functions and constants.
- /src/structs.rs - Data structures
- /src/data.rs - Process the input data
- /src/solver.rs - The actual search algorithm and evaluation function
- /src/bin/fireball.rs - The front-end (the actual program)
- /tests/integration_test.rs - Random data generator

## Data Format
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
