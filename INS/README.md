# INS Library

The MSS Inertial Navigation System (INS) library contains editable simulation examples and supporting functions for aided inertial navigation and attitude estimation. The implementations use high-rate inertial measurement unit (IMU) data with lower-rate aiding measurements and are compatible with MATLAB and GNU Octave.

## Simulation Examples

| File | Description |
| --- | --- |
| [`SIMaidedINSeuler.m`](SIMaidedINSeuler.m) | Compares Euler-angle error-state Kalman filter (ESKF) architectures using compass, position, optional velocity, or external attitude and heading reference system (AHRS) aiding. |
| [`SIMaidedINSquat.m`](SIMaidedINSquat.m) | Quaternion-based aided INS using a multiplicative error-state Kalman filter with position, optional velocity, and magnetometer or compass aiding. |
| [`SIMaidedINSheave.m`](SIMaidedINSheave.m) | Pressure-aided ESKF for estimating heave position, vertical velocity, and accelerometer bias. |
| [`SIMquatMEKF.m`](SIMquatMEKF.m) | Multiplicative extended Kalman filter (MEKF) example for quaternion attitude and angular-rate sensor bias estimation. |
| [`SIMquatObserver.m`](SIMquatObserver.m) | Nonlinear quaternion attitude observer with angular-rate sensor bias estimation. |

The simulation parameters, sensor frequencies, noise levels, aiding options, and test-signal selections can be edited near the beginning of each simulation file.

## Supporting Functions

The [`functions`](functions/) directory contains:

- `acc2rollpitch.m`: Static roll and pitch from accelerometer measurements, with optional bias compensation.
- `ins_euler.m`: 15-state compass- and position-aided INS ESKF using Euler angles.
- `ins_ahrs.m`: 9-state position-aided INS ESKF using attitude supplied by an external AHRS.
- `ins_mekf.m`: Quaternion-based INS ESKF aided by magnetometer and position measurements.
- `ins_mekf_psi.m`: Quaternion-based INS ESKF aided by compass and position measurements.
- `ins_heave.m`: Pressure-aided heave ESKF.
- `quatMEKF.m`: Quaternion MEKF for attitude and sensor-bias estimation.
- `quatObserver.m`: Nonlinear quaternion attitude observer.
- `insSignal.m`: Repeatable INS and IMU test-signal generator.
- `staticRollPitchYaw.m`: Static roll, pitch, and magnetic heading from accelerometer and magnetometer measurements.
- `gravity.m`: WGS-84 gravity model as a function of latitude.
- `magneticField.m`: Demonstration magnetic-field reference vectors and locations.

## Getting Started

Add MSS and its subfolders to the search path. In MATLAB, run the following command from the MSS root directory:

```matlab
mssPath
```

For GNU Octave, follow the [GNU Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20Octave.md).

Then run one of the simulation examples, for example:

```matlab
SIMaidedINSquat
```

## Reference

T. I. Fossen (2027). *Handbook of Marine Craft Hydrodynamics and Motion Control*, 3rd Edition, Wiley (in progress).
