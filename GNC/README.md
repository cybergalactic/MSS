# GNC Library

The GNC library contains reusable MATLAB and GNU Octave functions for marine guidance, navigation, and feedback control.

## Guidance and Path Following

- `LOSchi.m`, `ILOSpsi.m`, and `ALOSpsi.m`: two-dimensional line-of-sight guidance for course or heading control.
- `ALOS3D.m`: adaptive LOS guidance for three-dimensional heading and pitch control.
- `ALOSpsiHermite.m`: LOS guidance along a cubic Hermite-spline path.
- `LOSobserver.m`: filtered LOS-angle and LOS-rate estimation.
- `crosstrack.m`, `crosstrackWpt.m`, and `crosstrackWpt3D.m`: path-projection and tracking-error calculations.
- `hermiteSpline.m`, `hybridPath.m`, `getPathSignals.m`, and `projectToPath.m`: smooth path generation, evaluation, and projection.
- `addIntermediateWaypoints.m`, `order3.m`, and `order5.m`: waypoint preprocessing and polynomial path generation.

## Navigation and Attitude

- `EKF_5states.m`: estimates position, speed over ground, course over ground, and course rate from GNSS positions.

The more extensive inertial-navigation implementations are documented in the [INS Library](../INS/README.md).

## Control and Allocation

- `headingAutopilot.m`: SISO PID pole-placement heading controller with a
  third-order reference model.
- `PIDnonlinearMIMO.m`: nonlinear MIMO PID regulator for dynamic positioning.
- `integralSMCheading.m`: integral sliding-mode heading controller.
- `lqtracker.m`: linear-quadratic tracker design.
- `allocPseudoinverse.m`: unconstrained weighted control allocation.
- `refModel.m` and `refModelPolyExp.m`: command and reference-model generation.
- `sat.m` and `satlim.m`: symmetric and asymmetric signal saturation.
- `nomoto.m`: frequency-response plots for first- and second-order Nomoto models.

## Getting Started

Add MSS and its subfolders to the search path. In MATLAB, run the following command from the MSS root directory:

```matlab
mssPath
```

For GNU Octave, follow the [GNU Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20Octave.md).

The [`mssExamples`](../mssExamples/) and [`mssDemos`](../mssDemos/) directories contain examples that use these functions.
