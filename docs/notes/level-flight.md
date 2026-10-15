# Level flight

`spirula sfm auto --level-flight SLOPE`. Code: `src/sfm/map/LevelFlightPriors.h`, wired in
`src/sfm/Pipeline.cpp` (the mapper's prior source). Test: `sfm_level_flight_test`.

## What it is for

Straight-down photos over flat ground tie neighbouring flight lines together in height only
weakly. Without a position source (GPS, telemetry), a group of photos and the ground under it can
rise or sink together at almost no cost in reprojection, and each run settles somewhere different.
The ground then folds where the lines meet, and an orthophoto or DSM made from the model shows
steps that move from run to run.

Measured on a fixed-wing flight with no GPS (592 photos, rolling-shutter
correction on): reprojection about 1 px everywhere, no photo or image region worse than the rest,
but folds in the sparse ground of up to 5 m at a 50 m flight height. A photo's camera height moved
with the ground under it: 0.4 to 0.7 m per metre of fold, correlation 0.26 to 0.72 over five runs
that differed in one setting each. Turns, flight direction, the lens model, sequence pairs, feature
count and the rolling-shutter windows explained none of it.

## Model

An aircraft that holds its altitude puts cameras that are near each other over the ground at
nearly the same height. Per solve, from the posed cameras:

- up is against the cameras' mean viewing direction; when that is shorter than `min_nadir` (0.8)
  of the camera count, the cameras do not look one way and nothing is stated;
- each camera is tied to its `neighbours` (6) nearest over the ground **from other passes**: more
  than `min_gap` (10) photos away in the shooting order (the natural order of the file names), and
  to the side of the camera's own track (at least 30 degrees off its direction of travel, taken from
  the shorter of its steps to the photos before and after it, so a line's end does not take the
  jump to the next line as its direction). Along a line the photos are a few metres apart and the
  images fix their heights best; tying the nearest cameras regardless of pass tied each camera to
  its own line only, which fought the images and made the folds worse (measured, first version);
- one factor per pair: `up . (c_i - c_j) ~ 0`, sigma `SLOPE * distance over the ground`, Huber as
  every other prior. The level axes carry no weight.

The factors are relative: they state no position and leave the gauge alone. Another source's
factors (GPS, telemetry) pass through unchanged; the final gauge is fixed from that source, not
from these.

## Choosing SLOPE

`0.02` lets two cameras 50 m apart differ by 1 m at one sigma; past about 1.3 sigma the Huber loss
makes a real climb cost little more. A multirotor that follows the terrain, or a capture that is
not a survey, should not use it.
