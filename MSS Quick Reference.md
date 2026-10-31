# MSS Quick Reference

The m-files of the MSS (Marine Systems Simulator) toolbox are compatible with MATLAB (www.mathworks.com) and the free software GNU Octave (www.octave.org), facilitating broad accessibility and application in marine systems simulation. To use the MSS toolbox, please ensure the files are downloaded and correctly set up in your MATLAB/Octave environment. See the installation guides for [MATLAB](How%20to%20install%20MSS%20for%20MATLAB.md) and [GNU Octave](How%20to%20install%20MSS%20for%20GNU%20Octave.md).

This reference lists the user-facing functions, models, examples, demos, and
workflows. Internal callbacks, raw model-data scripts, and helper functions
used only by another documented entry point are not listed separately.

```matlab
>> mssHelp        % MSS Quick Reference 
>> mssPath        % Build/update and save the MSS path (not available in GNU Octave)
```

## Table of Contents
- [CRAFT (m-files)](#craft-m-files)
  - [Craft simulators](#craft-simulators)
  - [Craft models](#craft-models)
- [GNC (m-files)](#gnc-m-files)
- [LIBRARY (m-files)](#library-m-files)
  - [Modeling](#modeling)
  - [Kinematics](#kinematics)
  - [Environment](#environment)
  - [Ship maneuvers and visualization](#ship-maneuvers-and-visualization)
  - [Motion sickness](#motion-sickness)
  - [Transformations](#transformations)
  - [Numerical methods](#numerical-methods)
  - [Signal filters](#signal-filters)
- [INS (m-files)](#ins-m-files)
  - [INS error-state Kalman filter (ESKF) simulators](#ins-error-state-kalman-filter-eskf-simulators)
  - [Functions](#functions)
- [MSS Demos (m-files)](#mss-demos-m-files)
- [MSS Examples (m-files)](#mss-examples-m-files)
- [HYDRO](#hydro)
  - [Processing of data from hydrodynamic codes (m-files)](#processing-of-data-from-hydrodynamic-codes-m-files)
  - [Capytaine integration](#capytaine-integration)
  - [Vessel data (mat-files for hydrodynamic simulation)](#vessel-data-mat-files-for-hydrodynamic-simulation)
  - [Hydrodynamics (m-files)](#hydrodynamics-m-files)
  - [Hydrodynamic utilities (m-files)](#hydrodynamic-utilities-m-files)
  - [WAMIT geometry utilities (m-files)](#wamit-geometry-utilities-m-files)
- [Frequency-domain identification (FDI) of radiation models (m-files)](#frequency-domain-identification-fdi-of-radiation-models-m-files)
  - [FDI demos](#fdi-demos)
  - [Utils](#utils)
---

## CRAFT (m-files)

Start with the [CRAFT documentation](CRAFT/CRAFT-Library.md) for an overview of the
established craft models, their simulation scripts, and the WAMIT and
MSS-Capytaine workflow for user-defined vessels.

CRAFT includes established, ready-to-run mathematical models and a unified
workflow for vessels built from seakeeping data. `SIMhydroVessel` accepts the
common MSS `vessel` structure generated from WAMIT or Capytaine data.
WAMIT requires a commercial license; Capytaine is a free, open-source
alternative. Example vessel data and processing instructions are provided in
the [HYDRO](#hydro) catalog.

### Craft simulators
.../MSS/CRAFT/

```matlab
SIMclarke83     % Simulate a generic ship model characterized by L, B, and T under PD control.
SIMdsrv         % Simulate DSRV.m with an autopilot for depth control using successive-loop closure.
SIMfrigate      % Simulate frigate.m using a PID heading autopilot
SIMhydroVessel  % Simulate a 6-DOF vessel model generated from WAMIT or Capytaine data
SIMlauvMarie    % Simulate LAUV Marie with currents, directional waves, depth/heading control, and 3-D ALOS guidance
SIMmariner      % Simulate mariner.m with heading control and 2-D LOS straight-line path course control
SIMotter        % Otter heading/path-following simulation with currents, waves, GUI setup, and a legacy-model option
SIMosv          % Simulate osv.m under nonlinear DP control with constrained control allocation (dynamic optimization)
SIMcontainer    % Simulate container.m and Lcontainer.m under PD control
SIMnavalvessel  % Simulate navalvessel.m under PD control
SIMnpsauv       % Simulate npsauv.m with MIMO PID autopilots for depth and heading control, and 3-D straight-line path following using ALOS
SIMremus100     % Simulate remus100.m with autopilots for depth and heading control, and 3-D straight-line path following using ALOS
SIMsemisub      % Simulate the 6-DOF semisubmersible model under PID control
SIMsupply       % Simulate the linear supply vessel model under DP control
SIMtanker       % Simulate tanker.m under PD control
SIMzeefakkel    % Simulate zeefakkel.m using a PID heading autopilot
```

### Craft models
.../MSS/CRAFT/

```matlab
clarke83        % Linear maneuvering model parametrized using (L,B,T) found from linear regression of model tests (Clarke et al. 1983)
container       % Nonlinear maneuvering model of a high-speed container ship, L = 175 m, including roll (Son and Nomoto 1982)
DSRV            % Deep Submergence Rescue Vehicle (DSRV), L = 5.0 m (Healey 1992)
frigate         % Nonlinear autopilot model for a frigate, L = 100 m
hydroVessel     % Unified nonlinear 6-DOF model for WAMIT and Capytaine vessel structures
lauvMarie       % Unified nonlinear seakeeping/maneuvering model of LAUV Marie based on MSS-Capytaine data
Lcontainer      % Linearized model of a high-speed container ship, L = 175 m, including the roll mode (Son and Nomoto 1982)
mariner         % Nonlinear maneuvering model for the Mariner class vessel, L = 160 m 
navalvessel     % Nonlinear maneuvering model of a multipurpose naval vessel, L = 51.5 m
npsauv          % Naval Postgraduate School Autonomous Underwater Vehicle (AUV), L = 5.3 m (Healey and Lienard 1993)
osv             % Nonlinear model of an Offshore Supply Vessel (OSV), L = 83.0 m
otter           % Unified 6-DOF Otter USV maneuvering/seakeeping model using MSS-Capytaine data, L = 2.0 m
otter_legacy    % First-principles Otter USV model developed in 2019-2026
remus100        % REMUS 100 Autonomous Underwater Vehicle (AUV), L = 1.9 m
rig             % Semisubmersible linear mass-damper-spring model, L = 84.6 m
supply          % Linear DP model of a supply vessel, L = 76.2 m
tanker          % Nonlinear course unstable maneuvering model of a tanker, L = 304 m 
zeefakkel       % Nonlinear autopilot model of a recreational craft, L = 45 m
```

## GNC (m-files)
.../MSS/GNC/

```matlab
addIntermediateWaypoints % Add points along waypoint line segments for increased resolution
allocPseudoinverse   % Unconstrained control allocation
ALOS3D               % ALOS guidance laws for heading and pitch control in 3-D
ALOSpsi              % ALOS guidance law for heading control in 2-D
ALOSpsiHermite       % Computes the cross-track error and LOS angle to a cubic Hermite spline defined by waypoints
crosstrack           % Computes the path-tangential origin and cross-track error for a target
crosstrackWpt        % Computes the cross-track error when the path is a straight line between two waypoints
crosstrackWpt3D      % Computes the 3-D tracking errors (along-, cross- and vertical-track errors)
EKF_5states          % Estimation of SOG, COG, and course rate from NED positions or latitude-longitude
getPathSignals       % Evaluate a hybrid parameterized path and its derivatives
hermiteSpline        % Computes a cubic Hermite spline and the tangents to the spline for a given waypoint
hybridPath           % Generates coefficients for sub-paths between waypoints
integralSMCheading   % Integral sliding mode controller for heading control    
LOSchi               % LOS guidance law for course control in 2-D
LOSobserver          % Estimates the desired LOS angle and LOS rate from a LOS guidance law command
ILOSpsi              % ILOS guidance law for heading control in 2-D
lqtracker            % Computes the LQ tracker gain matrices for LTI systems
nomoto               % Generates Bode plots for the 1st- and 2nd-order Nomoto models
order3               % Path generation using cubic polynomials
order5               % Path generation using 5th-order polynomials
PIDnonlinearMIMO     % Nonlinear MIMO PID regulator for dynamic positioning (DP)
projectToPath        % Computes the orthogonal projection of a vehicle position onto a sampled 2-D path.
refModel             % Third-order reference model for position, velocity, and acceleration
refModelPolyExp      % Transition from x_start to x_final using a polynomial or an exponential curve
sat                  % Symmetric elementwise signal saturation
satlim               % Asymmetric elementwise signal saturation
controlMethod        % GUI for selecting a control method from a list
displayVesselStructure % Print the fields and values of an MSS vessel structure
isoctave             % Return true when running under GNU Octave
rangeCheck           % Validate that values lie within inclusive limits
```


## LIBRARY (m-files)
.../MSS/LIBRARY/

### Modeling

```matlab
addedMassSurge      % Hydrodynamic added mass in surge, A11, approximated by the formula of Söding (1982)
coeffLiftDrag       % Hydrodynamic lift and drag coefficients as a function of angle of attack of a submerged "wing profile"
forceLiftDrag       % Hydrodynamic lift and drag forces as a function of the angle of attack of a submerged "wing profile" 
crossFlowDrag 	    % Crossflow drag computed from strip theory integrals
Dmtrx               % 6x6 linear damping matrix for marine craft (submerged and floating)
Gmtrx               % 6x6 system spring stiffness matrix G
GM_surfaced2submerged % GM_T and BM_T transition for an AUV diving from the surface to a given depth
gvect               % 6x1 vector of restoring forces, Euler angles as input
gRvect              % 6x1 vector of restoring forces, Euler angle, or unit quaternion rotation matrix as input
imlay61             % 6x6 hydrodynamic added mass and Coriolis-centripetal matrices MA and CA for a prolate spheroid
m2c                 % 6x6 Coriolis-centripetal matrix C(nu) from system inertia matrix M
rbody               % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a general body
spheroid            % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a prolate spheroid 
thrConfig           % 3xr thruster configuration matrix for main propellers, tunnel thrusters, and azimuth thrusters
vesselPeriods       % Periods, natural frequencies, and damping ratios from frequency-dependent vessel data
wageningen          % Thrust and torque coefficients of the Wageningen B-series propellers 
```

### Kinematics

```matlab
ecef2llh            % Longitude, latitude, and height from ECEEF positions x, y, and z
euler2q             % Unit quaternion from Euler angles
eulerang            % Computes the Euler angle transformation matrices J, Rzyx and Tzyx
flat2llh            % Longitude, latitude, and height from flat-earth positions x, y, and z
llh2ecef            % ECEEF positions x, y, and z from longitude, latitude, and height
llh2flat            % Flat-Earth positions x, y, and z from longitude, latitude, and height
modalFromMDG        % Compute modal properties (time const., natural freqs., damping ratios) of a linear 6-DOF model
R2euler             % Euler angles from rotation matrix elements
Rll                 % Euler angle rotation matrix Rll for longitude and latitude
Rquat               % Unit quaternion rotation matrix R in SO(3)
Rzyx                % Euler angle rotation matrix R in SO(3)
Tquat               % Unit quaternion transformation matrix T, representing the attitude dynamics
Tzyx                % Euler angle transformation matrix T, representing the attitude dynamics
q2euler             % Euler angles from a unit quaternion
quatern             % Unit quaternion transformation matrix J
quatprod            % Quaternion product
quest               % Quaternion rotation matrix R(q) and unit quaternion q between two vectors W = R(q) V
quest6dof           % 6-DOF vector eta = [x,y,z,phi,theta,psi] from three marker positions using the QUEST algorithm
```

### Environment

```matlab
blendermann94       % Compute the wind forces and wind coefficients using Blendermann (1994)
encounter           % Encounter frequency as a function of wave peak frequency, vessel speed, and wave direction
hs2vw               % Convert significant wave height into an equivalent wind speed
isherwood72         % Compute the wind forces and coefficients based on Isherwood (1972) 
spheroidRAO         % Compute 6-DOF RAOs for a submerged prolate-spheroid-shaped AUV in waves.
torsetSpectrum      % Torsethaugen double-peaked wave spectrum
vw2hs               % Converts average wind speed to significant wave height
waveEncounterPhase  % Accumulate signed encounter phase during maneuvering using explicit waveModel history
waveForceRAO        % Compute 1st-order wave loads and elevation; return updated waveModel, with optional depth attenuation
waveInitialization  % Initialize the spectrum, random phases, RAO tables, encounter-phase history, and reference depth
waveMotionRAO       % Compute WF position, velocity, acceleration, and elevation; return updated waveModel
waveresponse345     % Steady-state heave, roll, and pitch responses for a ship in regular waves 
waveSpectrum        % Function computing state-of-the-art wave spectra
waveDirectionalSpectrum % Computes the directional wave spectrum using a spreading function
```

### Ship maneuvers and visualization

```matlab
animateShip         % Animates a Viking ship moving along a specified path on a North-East (y-x) plot
pullout             % Ship pullout maneuver
turncircle          % Ship turning circle
zigzag              % Zigzag maneuver for 3-DOF models
zigzag6dof          % Zigzag maneuver for 6-DOF models
```

### Motion sickness

```matlab
ISOmsi              % ISO 2631-3, 1997 motion sickness incidence
HMmsi               % O'Hanlon and McCauley (1974) motion sickness incidence
```

### Transformations

```matlab
conversion          % Defines global conversion factors for GNC applications
ssa                 % Smallest signed angle maps an angle in rad to the interval [-pi pi) or [-180 180)
Smtrx               % 3x3 vector skew-symmetric matrix S
Hmtrx               % 6x6 system transformation matrix H
vex                 % Computes a = vex(S(a)) where S is a skew-symmetric matrix
```

### Numerical methods

```matlab
euler2              % Integrates a system of ordinary differential equations using Euler’s 2nd-order method
expm_squaresPade    % Custom-made matrix exponential using Pade approximations (exponential map for matrix Lie groups)
expm_taylor         % Custom-made matrix exponential using Taylor series (exponential map for matrix Lie groups)
invQR               % Matrix inversion using QR factorization for improved numerical stability
rk4                 % Integrates a system of ordinary differential equations using Runge-Kutta’s 4th-order method
```

### Signal filters

```matlab
highPassFilter      % First-order high-pass filter using exact discretization
lowPassFilter       % First-order low-pass filter using exact discretization
makeIntegrator      % Create a stateful discrete-time integrator function handle
makeLowPass         % Create a stateful first-order low-pass filter function handle
notchFilter         % Second-order notch filter using RK4 discretization or IIR filtering
sawToothWave        % Sawtooth wave signal, which can be used for testing
waveFreqObserver    % Wave encounter frequency estimator (Belleter, Galeazzi and Fossen 2015)
```

## INS (m-files)
.../MSS/INS/

### INS error-state Kalman filter (ESKF) simulators

```matlab
SIMaidedINSeuler     % Simulate the ESKF for aided INS using Euler angles
SIMaidedINSheave     % Simulate the ESKF for aided INS in heave using pressure measurements
SIMaidedINSquat      % Simulate the ESKF for aided INS using unit quaternions (MEKF representation)
SIMquatMEKF          % Simulate the MEKF quaternion-based attitude observer for 9-DOF IMU measurements
SIMquatObserver      % Simulate the nonlinear quaternion-based attitude observer for 9-DOF IMU measurements
```

### Functions

```matlab
acc2rollpitch        % Static roll and pitch angles from IMU specific force measurements
gravity              % Acceleration of gravity as a function of latitude using the WGS-84 ellipsoid parameters
ins_ahrs             % Error-state Kalman filter (ESKF) for an INS aided by position and AHRS measurements 
ins_euler            % Error-state Kalman filter (ESKF) for an INS aided by position and yaw angle measurements
ins_heave            % Error-state Kalman filter (ESKF) for an INS (heave only) aided by pressure measurements
ins_mekf             % Error-state Kalman filter (ESKF) for an INS aided by position and magnetic field measurements
ins_mekf_psi         % Error-state Kalman filter (ESKF) for an INS aided by position and yaw angle measurements
insSignal            % INS signal generator for testing of Kalman filters and observers
magneticField        % NED magnetic field reference vectors 'm_ref' for different cities
quatObserver         % Nonlinear quaternion-based attitude observer for 9-DOF IMU measurements
quatMEKF             % MEKF quaternion-based attitude observer for 9-DOF IMU measurements
staticRollPitchYaw   % Static roll, pitch, and magnetic heading from IMU specific force and magnetometer measurements
```

## MSS demos (m-files)
/MSS/mssDemos/

```matlab
1) KinDemo     % Euler angle and quaternion kinematics
2) ManDemo     % Maneuvering trials
3) StabDemo    % Straight-line, directional, and positional motion stability
4) WaveDemo    % Wave spectra demonstrations
```

## MSS examples (m-files)
/MSS/mssExamples/

```matlab
exAUVHydrostatics  % Computation of the hydrostatic quantities for a cylinder-shaped AUV
exBoxShapedShip    % Computation of the transverse metacentric height and the heave/roll periods of a box-shaped ship
exEKF              % Discrete-time extended Kalman filter for a nonlinear second-order system
exEKFparamterEstimation % Extended Kalman filter for simultaneous state and parameter estimation
exFeedback         % For-loop implementation for numerical integration of a 1st-order system under feedback and feedforward control
exFFT              % Estimation of the wave encounter frequency from time series using the fast-Fourier transform (FFT)
exPlotGM           % Compute and plot the GM_T and BM_T for an AUV diving from the surface to a given depth
exHermiteWaypoints % Computes a C1 Hermite/MAKIMA spline path through waypoints given in North-East coordinates.
exHybridPath       % Computation of a hybrid continuous path parametrized by waypoints
exINS_AHRS         % Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and AHRS attitude measurements 
exINS_Euler        % Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and compass measurements
exINS_MEKF         % Unit quaternion error-state (indirect) Kalman filter for INS aided by position and magnetic field measurements 
exINSwaveFilter    % Wave filtering techniques for Inertial Navigation System (INS) measurements
exIntWindup        % Demonstrates integrator windup and anti-windup when the control law is saturated.
exKF               % For-loop implementation (predictor-corrector representation) of a discrete-time linear Kalman filter (KF) 
exKalmanWavefilterAutopilot % Discrete-time Kalman wave filtering and PID heading control of the Mariner
exKT               % Computation of the Nomoto gain K and time constant T from a step response using nonlinear least-squares
exLinspec          % Linear approximations to the PM, JONSWAP, and Torsethaugen spectra using nonlinear least-squares
exLQFinHor         % LQ finite time-horizon tracking controller for a mass-damper-spring system
exLQR              % Computes the LQR gains for a mass-damper system
exLQtrack          % Computes the LQ optimal tracking gains for a mass-damper system
exManeuveringModel % Maneuvering model from frequency-dependent hydrodynamic coefficients A(ω) and B(ω), and the wave spectrum S(ω).
exMDS              % Plots the step response of a 2nd-order mass-damper system
exMSI              % Plots the ISO 2631-1 (1997) and O'Hanlon and McCauley (1974) Motion Sickness Incidence (MSI) curves
exNomoto           % Bode plots of ships parametrized by Nomoto’s 1st- and 2nd-order models
exObsCtr           % Observability and controllability matrices of a supply vessel
exOtter            % EKF/course-control example using the first-principles otter_legacy model
exPassiveObs       % Plots the loop transfer function of the passive observer used for heading control
exPassiveWavefilterAutopilot % Passive wave filtering and PID heading control using estimated or measured yaw rate
exPassiveWavefilterDP % Full-scale passive DP wave filtering with GNSS, compass and nonlinear MIMO PID control
exPathGen          % Path generation using cubic polynomials 
exPlotRAO          % Script for plotting motion and force RAOs
exPullout          % Performs a pullout maneuver for two different ships
exQuadProg         % Quadratic programming applied to waypoint trajectory generation
exQuest            % 6-DOF position/attitude vector from camera measurements using the QUEST algorithm
exRefMod           % 2nd-order reference model with nonlinear damping and velocity saturation
exLoadCondition    % Plot heave, roll, and pitch periods versus draft and metacentric height
exRegularWaveEstimation % Compare manual heave, roll and pitch excitation estimates with MSS-Capytaine force RAOs for supply and Otter (80 kg)
exResonance        % Computes the closed-form responses in heave, roll, and pitch for a marine craft exposed to regular waves
exRRD1             % Roll and sway-yaw transfer functions for the Son and Nomoto container ship
exRRD2             % Rudder-roll damping (RRD) system for the Son and Nomoto container ship
exRRD3             % Inverse response in roll for the Son and Nomoto container ship due to a right-half-plane zero  (non-minimum phase)  
exShipHydrostatics % Computation of the hydrostatic quantities for a ship 
exSMC              % Integral sliding mode control (SMC) design for heading control
exSpline           % Path generation using cubic Hermite spline interpolation 
exSTA              % Adaptive-gain super twisting algorithm (STA) for heading control
exTurnCircle       % Generates the turning circle for two different ships
exWageningen       % Computes thrust and torque curves for a propeller using the Wageningen B-series data
exWaveFreqObserver % Nonlinear observer for estimation of the wave encounter frequency
exWaveForceRAO     % Wave elevation and generalized 1st-order wave forces from force RAOs using different wave spectra 
exWaveMotionRAO    % Wave elevation and ship wave-frequency (WF) motions from motion RAOs using different wave spectra
exWindForce        % Plots the wind coefficients by Isherwood (1972) 
exWOPC             % Fully and underactuated weather-optimal positioning control
exZigZag           % Generates zigzag maneuvers for two different ships and the Remus 100 AUV
```

## HYDRO
/MSS/HYDRO/

### Processing of data from hydrodynamic codes (m-files)

```matlab
vessel2ss          % computes the fluid-memory transfer functions and stores the data as a mat-file containing the structure <vesselABC>
wamit2vessel       % Reads data from WAMIT output files and stores the data as a mat-file containing the structure <vessel>
```

The WAMIT importer and MSS-Capytaine exporter create the frequency-dependent vessel
structure. The power-based maneuvering model is then computed separately.
Floating vessels use three target time constants in seconds and three
damping-ratio increments, while submerged vehicles use four target time
constants in seconds and two damping-ratio increments:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,plotFlag);
vessel = computeManeuveringModel(vessel,omega_p,T_1236,delta_zeta_45,plotFlag);
```

`T_126 = [T_surge T_sway T_yaw]` specifies total diagonal velocity-decay
targets, with `Bv(i,i) = M(i,i)/T_i - B_eq(i,i)`. Values must be positive
and must not require negative additional viscous damping. The surface-vessel
GUI displays these time constants in seconds. Stored MAT-file inputs supply
the defaults for all catalogue vessels. `wamit2vessel` collects and saves
these inputs, with Return accepting its vessel-specific proposals.
Averaging uses JONSWAP with gamma 3.3 and a fixed peak `omega_p = 1.0` rad/s
when the peak argument is omitted or empty. Pass another peak explicitly.

Old saved `vessel.powerBased.kappa_126` percentages are reused with a warning
only when the damping arguments are omitted. Explicit three-element damping
arguments always mean `T_126` in seconds; do not reuse percentage vectors.

This adds `vessel.powerBased`, including the equivalent matrices `A_eq` and
`B_eq` and the viscous damping matrix `Bv`, together with the maneuvering-model
matrices `M`, `D`, and `G`. Use `vessel.powerBased.Bv` for additive viscous damping; the former
top-level field `vessel.Bv` is not created.

See the processing instructions for the [WAMIT data](HYDRO/vessels_wamit/WAMIT-Vessels.md)
and [Capytaine data](HYDRO/vessels_capytaine/Capytaine-Vessels.md).

### Capytaine integration

[MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine) is a separate Python add-on that uses the open-source [Capytaine](https://capytaine.org/) boundary-element solver and exports the results in the MSS `vessel` structure format. Clone MSS-Capytaine separately and run `python main.py` from its project directory to regenerate the example hydrodynamic data.

Pre-generated [testShip](HYDRO/vessels_capytaine/testShip/Test-Ship.md) and submerged
[LAUV_marie](HYDRO/vessels_capytaine/LAUV_marie/LAUV-Marie.md) cases are included in MSS.
Python and Capytaine are required only to regenerate the data; the included
mat-files can be loaded and processed directly in MATLAB or GNU Octave.
The dedicated [`lauvMarie.m`](CRAFT/AUV/models/lauvMarie.m) and
[`SIMlauvMarie.m`](CRAFT/AUV/SIMlauvMarie.m) pair uses the LAUV data as a
standalone, wave-enabled craft model with vehicle-specific nonlinear and
calibratable terms. `SIMhydroVessel` provides a generic hydrodynamic test tool:
use **Load vessel MAT-file** to load exports from
`MSS-Capytaine/vessels_capytaine/<vessel>_<format>/`. Test ship, Supply and S175
are bundled Capytaine cases; Tanker, FPSO and Semisubmersible are WAMIT cases. LAUV Marie can be
loaded manually; use `SIMlauvMarie` for its dedicated simulation.
Their JSON-selected linear damping inputs are stored in `vessel.powerBased`
and used as the defaults by `hydroVesselConfig.m`. Simulation-specific overrides can be set there. All three WAMIT vessels also store their selected damping inputs.

### Vessel data (mat-files for hydrodynamic simulation)

```matlab
testShip             % Capytaine data for the synthetic MSS-Capytaine test ship
LAUV_marie           % Capytaine data for NTNU's submerged LAUV Marie case
fpso, fpsoABC       % WAMIT data for a FPSO
semisub, semisubABC % WAMIT data for a semisubmersible
tanker, tankerABC   % WAMIT data for a tanker
s175                 % Capytaine zero-speed S175 container ship
supply               % Capytaine zero-speed supply vessel
otter_55kg, otter_80kg % Capytaine Otter USV load conditions
wigley               % Capytaine benchmark hull
s175_speedDependent  % Legacy S175, only in mssExamples/data for exManeuveringModel
```

After loading the data files to the workspace using the Matlab command load, the following data structures are available:

| vessel   | vessel.main | vesselABC |
| -------- | --------    | --------  |
|          main: [1×1 struct]      |      name: 'tanker' | Ar: {6×6 cell} |
|    velocities: [1×N_speed]              |         T: draft    | Br: {6×6 cell} |
|      headings: [1×36 double]     |         B: beam     | Cr: {6×6 cell} |
|           MRB: [6×6 double]      |       Lpp: length   | Dr: {6×6 cell} | 
|             C: [6×6×60 double]   |         m: mass     | MRB: [6×6 double] | 
|         freqs: [1×60 double]     | rho: density of water | MA: [6×6 double] | 
|             A: [6×6×60 double]   | k44: radius of gyration | G: [6×6 double] | 
|             B: [6×6×60 double]   | k55: radius of gyration | Minv: [6×6 double] | 
|     motionRAO: [1×1 struct]      | k66: radius of gyration | r_g: [x_g y_g z_g] | 
|      forceRAO: [1×1 struct]      | g: acceleration of gravity | |
|         drift: [1×1 struct]      | nabla: volume displacement | |
|    powerBased: [1×1 struct]      | CB: center of buoyancy | |
| roll: [1×1 struct] (legacy S175)       | GM_T: transverse metacentric height | |
|                                  | GM_L: longitudinal metacentric height | |
|                                  | CG: center of gravity | |

All Capytaine and WAMIT catalogue files contain selected damping inputs
under `vessel.powerBased`; `computeManeuveringModel` adds or recomputes the
equivalent matrices and `Bv`. The retained legacy S175 also contains its
original `vessel.roll` structure with Ikeda damping `vessel.roll.Bv44`.
Mean drift uses only `vessel.drift.w` and `vessel.drift.coeff{1:3}` for
surge/sway/yaw, in loads per squared incident-wave amplitude. WAMIT moments
retain their source body origin; Capytaine exports use CG. `plotWD` does not
read the removed `vessel.driftfrc` fields.

The ShipX importer and catalogue have been removed. Supply and S175 now use
Capytaine data. The speed-dependent S175 required by the 2025 comparison is
retained only in [`mssExamples/data`](mssExamples/data/README.md).
`exManeuveringModel` uses a common JONSWAP baseline, identical additive
viscous damping in both methods, and a 300 s default duration.

Run `mssCraft` to display both the hydrodynamic MAT-file catalogue and all
CRAFT dynamics models, with companion simulation scripts and locations.


### Hydrodynamics (m-files)

```matlab
computeManeuveringModel % Computes power-based A_eq, B_eq, and Bv from A(ω) and B(ω)
cylinderDrag        % Reynolds- and aspect-ratio-dependent 2-D cylinder drag coefficient
Hoerner             % 2-D Hoerner cross-flow form coefficient as a function of B and T
plotAB_eq           % Plots equivalent Aij and Bij values (maneuvering model approximation)
plotABC             % Plots the hydrodynamic coefficients Aij, Bij, and Cij as a function of frequency 
plotBv              % Plots potential, equivalent, and power-based viscous damping
plotTF              % Plots the motion or force RAO transfer functions
plotWD              % Plots the wave drift amplitudes
```

### Hydrodynamic utilities (m-files)

```matlab
ABCtransform        % Transform and complete ShipX hydrodynamic matrices in MSS axes
imp2svd             % Singular values of the Hankel matrix of a sampled impulse response
retardation         % Compute a radiation retardation function from damping data
retardation2ss      % Estimate a state-space realization from a retardation function
```

### WAMIT geometry utilities (m-files)

```matlab
gdf2body            % Transform a WAMIT GDF geometry to the body-fixed origin
gdf_edit            % Repair an open waterline in a low-order WAMIT GDF file
gdf_waterline       % Remove panels above the waterline from a WAMIT GDF file
plot_tecplot        % Plot a WAMIT-generated Tecplot panel geometry
plot_wamitgdf       % Plot a WAMIT low-order GDF geometry
plot_wamitidf       % Plot a WAMIT intermediate-order IDF geometry
```

## Frequency-domain identification (FDI) of radiation models (m-files)

### FDI demos
```matlab
Demo_FDIRadMod_NA     % FDI using hydrodynamic data without infinite-frequency added mass
Demo_FDIRadMod_WA     % FDI using hydrodynamic data, including infinite-frequency added mass
```

### Utils
```matlab
EditAB                 % Function: preparing the data for identification, selecting frequency range, and removing wild-points
FDIRadMod              % Identify the SISO transfer function corresponding to the coupling specified
fit_siso_fresp         % Fit a continuous SISO transfer function to the frequency response data
ident_retardation_FD   % Identification of a parametric radiation convolution model K(s) = P(s)/Q(s)
ident_retardation_FDna % Identification of a parametric model A(jw) = B(w)/(jw) - A(w) 
```
