# MSS General Library

The MSS general library contains reusable MATLAB and GNU Octave functions for marine-craft modeling, kinematics, environmental loads, maneuvering tests, numerical integration, signal processing, and motion-sickness analysis.

## Catalogue

| Directory | Contents |
| --- | --- |
| [`modeling`](modeling/) | Rigid-body and hydrodynamic matrices, restoring forces, lift and drag, thruster configurations, propeller data, modal analysis, and vessel periods. |
| [`kinematics`](kinematics/) | Euler-angle and quaternion transformations, rotation matrices, skew-symmetric matrices, geodetic transformations, and attitude determination. |
| [`environment`](environment/) | Wind loads, encounter frequency, wave spectra, directional spreading, wave initialization, force and motion RAOs, and regular-wave response. |
| [`maneuvering`](maneuvering/) | Turning-circle, zigzag, pullout, and ship-animation utilities. |
| [`numericalMethods`](numericalMethods/) | Euler and Runge–Kutta integration, matrix-exponential approximations, and QR-based matrix inversion. |
| [`signalFilters`](signalFilters/) | Low-pass, high-pass, notch, wave-frequency, integration, and waveform utilities. |
| [`motionSickness`](motionSickness/) | ISO and O’Hanlon–McCauley motion-sickness-incidence calculations. |

## Selected Entry Points

- `rbody.m`, `m2c.m`, `Dmtrx.m`, and `Gmtrx.m`: standard 6-DOF system matrices.
- `crossFlowDrag.m`, `forceLiftDrag.m`, and `XuuITTC.m`: nonlinear hydrodynamic loads.
- `eulerang.m`, `Rzyx.m`, `Rquat.m`, and `Tquat.m`: attitude kinematics.
- `waveSpectrum.m`, `waveDirectionalSpectrum.m`, `waveInitialization.m`, `waveForceRAO.m`, and `waveMotionRAO.m`: irregular-wave modeling and response.
- `turncircle.m`, `zigzag.m`, and `zigzag6dof.m`: standard maneuvering trials.
- `rk4.m`: fourth-order Runge–Kutta integration.
- `lowPassFilter.m`, `highPassFilter.m`, and `waveFreqObserver.m`: signal filtering and wave-frequency estimation.

## Getting Started

In MATLAB, add MSS and its subfolders to the path by running the following command from the MSS root directory:

```matlab
mssPath
```

For GNU Octave, follow the [GNU Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20Octave.md).

See the [MSS Quick Reference](../MSS%20Quick%20Reference.md) for a function-by-function catalogue and the individual function headers for syntax, assumptions, and references.
