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

- `mssConstants.m`: common seawater density, air density, and standard gravity values used throughout MSS.
- `rbody.m`, `m2c.m`, `Dmtrx.m`, and `Gmtrx.m`: standard 6-DOF system matrices.
- `crossFlowDrag.m`, `forceLiftDrag.m`, and `XuuITTC.m`: nonlinear hydrodynamic loads.
- `eulerang.m`, `Rzyx.m`, `Rquat.m`, and `Tquat.m`: attitude kinematics.
- `waveSpectrum.m`, `waveDirectionalSpectrum.m`, `waveInitialization.m`, `waveForceRAO.m`, and `waveMotionRAO.m`: irregular-wave modeling and response.
- `waveEncounterPhase.m`: signed encounter-phase accumulation shared by the force- and motion-RAO functions.
- `turncircle.m`, `zigzag.m`, and `zigzag6dof.m`: standard maneuvering trials.
- `rk4.m`: fourth-order Runge–Kutta integration.
- `lowPassFilter.m`, `highPassFilter.m`, and `waveFreqObserver.m`: signal filtering and wave-frequency estimation.

## Wave Simulation

Initialize the spectrum and RAO tables once, then retain the updated
`waveModel` at each simulation update:

```matlab
[environment, waveModel] = waveInitialization(vessel, environment);
for t = 0:0.1:100
    [tau_wave1, waveElevation, waveModel] = waveForceRAO( ...
        t, waveModel, U, psi, beta_wave);
end
```

For motion RAOs, select the motion model and retain the fifth output:

```matlab
[environment, waveModel] = waveInitialization(vessel, environment, 'motion');
for t = 0:0.1:100
    [eta_WF, nu_WF, nudot_WF, waveElevation, waveModel] = waveMotionRAO( ...
        t, waveModel, U, psi, beta_wave);
end
```

The signed encounter rate is integrated with the trapezoidal rule instead
of multiplying its instantaneous value by elapsed time. This preserves phase
continuity during changes in speed and heading, including negative encounter
rates when overtaking wave components in following seas. The nonnegative
encounter-frequency magnitude remains the absolute value of this rate, as
described in Chapter 10 of Fossen (2027). Motion velocities and accelerations
use the quasi-steady RAO derivatives at the current operating point.

`beta_wave` is the propagation direction in NED; the relative direction is
`beta_wave - psi`. Pass the speed used in the book's forward-speed encounter
relation. Calls must be in nondecreasing time order. Each realization retains
its own phase history in `waveModel`; reinitialize it to restart. No persistent
wave state or `clear` statement is required.

**Updating existing simulation loops:** the input arguments are unchanged.
Retain the third output of `waveForceRAO` or the fifth output of
`waveMotionRAO` to accumulate phase when operating conditions change. Calls
that omit this output and use a fresh initialized model remain suitable for
direct evaluation at constant speed and heading. `SIMotter`, `SIMlauvMarie`,
`SIMhydroVessel`, and the RAO examples retain the updated model explicitly.

`waveForceRAO.m` also accepts an optional positive-down vehicle depth. For a
submerged vessel, `waveInitialization.m` reads the RAO reference depth from
`vessel.main.submergenceDepth`, and each force component is attenuated below
that reference by `exp[-(omega^2/g) Delta_d]`. Omitting the depth preserves
the original surface-vessel behavior, and shallower depths do not amplify the
loads.

## Getting Started

In MATLAB, add MSS and its subfolders to the path by running the following command from the MSS root directory:

```matlab
mssPath
```

For GNU Octave, follow the [GNU Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20Octave.md).

See the [MSS Quick Reference](../MSS%20Quick%20Reference.md) for a function-by-function catalogue and the individual function headers for syntax, assumptions, and references.
