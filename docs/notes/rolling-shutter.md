# Rolling shutter from neighbour matches

`spirula sfm rs`, and `spirula sfm auto --rolling-shutter`. Code: `src/sfm/rs/RollingShutter.h`
(the estimate), `src/sfm/rs/Apply.h` (applying it to feature files and photos),
`src/app/cli/sfm_rs.cpp` (the command).

## What it is for

A camera with a rolling-shutter sensor on a moving aircraft exposes each image row at a slightly
different moment. If the aircraft rotates during the read-out, every row sees the ground through a
slightly different orientation, and the photo is no longer the image of any pinhole camera. Over
flat ground the mapper explains the distortion with whatever it has: a focal length several times
too long, cameras tilted far off vertical, fragments that will not merge.

Measured on a camera with 4056 x 3040 resolution and 87 ms read-out on two fixed-wing
aircraft: pixels move by a median 32 to 57 px, and up to 157 px, within one photo.

## Model

Per photo, the camera's rotation during read-out is a cubic in the image row:

    s      = (row - height / 2) / height                      in -0.5 .. 0.5
    w(s)   = w1 s + w2 (s^2 - 1/12) + w3 (s^3 - 0.15 s)        three axes each, radians, camera frame
    ray    = Exp(w(s)) K^-1 [u v 1]

`w1` is a constant rate; `w2` and `w3` let the rate change within the read-out, and are needed:
with `w1` alone neighbouring photos agree to 2 to 3 px, with all three to about 1 px (the match
noise). The two higher polynomials are orthogonal to the lower ones over the frame.

Rotation does not depend on scene depth, so the correction is a pure remapping of pixels. The
aircraft's own translation during read-out is ignored: under 1% of the effect at these speeds.

## Estimation

No reconstruction exists yet when this runs, so the estimate is made on a deliberately simple
scene: flat ground, cameras looking roughly down.

1. Sequence order is the natural order of the file names (`photo9` before `photo10`).
2. Unverified matches between photos up to `max_gap` apart are filtered with a loose homography
   (10 px). A tight gate would remove exactly the distortion being measured.
3. Where consecutive photos share fewer than `min_link` matches the sequence is cut. Each unbroken
   run is split into overlapping windows of 60 photos, solved independently and in parallel.
4. A window starts from neighbours chained by a rigid fit, every camera level at one height.
5. Levenberg-Marquardt over, per photo, a ground shift (2), height (1), attitude (3) and the nine
   read-out terms; per window, a camera mounting rotation (3) and one radial distortion term. The
   residual is the mismatch of a match on the ground plane, in pixels.
6. Every photo takes its read-out terms from the window in which it is furthest from an edge.

Everything has a prior, and their balance matters more than their values:

| Unknown | Prior | Why |
|---|---|---|
| read-out, constant rate | 6 mrad over the frame | over a plane, a shear and stretch are also explained by pose; left loose they trade against height |
| read-out, higher orders | 30 mrad | a pose cannot imitate them |
| attitude | 3 degrees from level, 8 in yaw | |
| height | 1.5 m per photo step | the ground is smooth |
| matches | 15 px, Huber at 4 px | their errors are correlated; at their 1 px noise they overpower every prior and a strip drifts |

The focal length is not refined (over a plane it trades against height). When not given it is
scanned over a grid on the three best-matched windows; the fit has a clear minimum. The estimate is
insensitive to it: a focal length wrong by 1.5 times changes the correction by about as much as two
independent estimates of it differ.

## Applying it

- **Keypoints.** `correctPoint` moves each keypoint to where an instant exposure at the mid-row
  pose would have put it. Lens distortion is left in, for the mapper to estimate as usual. The
  corrected feature directory carries `rolling_shutter.txt`, which lists the estimate per photo
  and marks the directory as corrected.
- **Photos.** `--write-images` (or `auto --rolling-shutter-images`) resamples the photos the same
  way, bicubic, and copies their EXIF. Training and any other tool must use these, not the
  originals: the poses describe the corrected geometry.

In `auto` the stage sits between extraction and matching: one unverified matching pass, the
estimate, the feature files rewritten in place, then the pipeline as it was.

## Measured

Source build at `b9c8042` plus this stage, RTX 5090. Plain `auto` against `auto --rolling-shutter`,
no flight log, default settings:

| Capture | Photos | Plain | With the stage |
|---|---|---|---|
| second aircraft, 8 mm lens | 40 | 31 placed, 1.31 px, 22 s | 40 placed, 0.63 px, 12 s |
| forest, 8 mm lens with 13% distortion, two pieces | 55 | 25 placed, 1.19 px, 37 s | 37 + 18 placed, 0.76 px, 18 s |
| narrow 16 mm lens, GPS in EXIF | 111 | 69 placed, 2.00 px, 32 s | 84 placed (+ a second model), 0.59 px, 38 s |
| second aircraft, whole flight with a dozen turns | 1,041 | 88 placed in the largest of 25 models, 1.48 px, 392 s | 1,037 placed in one model, 0.70 px, 324 s |

The stage itself takes 4 to 9 s on the small sets and 17 s on 1,041 photos (22 windows). On the narrow lens the focal length found goes from
54,444 px to 10,443 px (true: about 10,400).

## Limits

- Roughly straight-down photos over ground that is flat within one photo. Oblique views, low
  flights over tall objects and steep terrain are outside the model.
- A sequence: consecutive file names must be consecutive photos, all one size.
- Runs shorter than 6 photos are left uncorrected.
- Banked turns fit worse (2 to 5 px against 1 px on straight stretches): the level-camera prior is
  wrong there.
- Photos with an EXIF orientation other than 1 are not handled by `--write-images`.
- Only for rolling-shutter sensors. On a global-shutter camera the fit still finds a "correction"
  and hurts the model: on a 1,001-photo nadir flight from about 700 m, reprojection error went
  from 0.71 to 0.96 px and the strips bent against each other by tens of metres.
- It does not make the mapper merge models that disagree about tilt. A narrow lens over flat
  ground leaves tilt and position ambiguous; that needs position and attitude priors.
