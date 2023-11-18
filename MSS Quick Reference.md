# MSS Quick Reference

A quick reference guide for the MATLAB MSS toolbox:

```matlab
>> help MSS               % List mss commands
>> help mssExamples       % List book examples (T. I. Fossen. Handbook of Marine Craft Hydrodynamics and Motion Control, Wiley 2021)
>> mssSimulink            % Simulink Library
```

## Table of Contents
- [Simulink demos](#simulink-demos)
- [GNC demos (m-files)](#GNC-demos-m-files)
- [Examples (m-files)](#examples-m-files)
- [Marine craft simulator (m-files)](#marine-craft-simulator-m-files)
  - [Marine craft models](#marine-craft-models)
  - [Time-series simulation](#time-series-simulation)
- [Utils (m-files)](#utils)
  - [Modelling](#modelling)
  - [Kinematics](#kinematics)
  - [Environment](#environment)
  - [Ship maneuvers](#ship-maneuvers)
  - [Motion sickness](#motion-sickness)
  - [Transformations](#transformations)
  - [Numerical integration methods](#numerical-integration-methods)
- [GNC (m-files)](#gnc-m-files)
  - [Guidance](#guidance)
  - [Navigation](#navigation)
  - [Control](#control)  
- [HYDRO](#hydro)
  - [Hydrodynamic templates (Simulink)](#hydrodynamic-templates-simulink)
  - [Processing of data from hydrodynamic codes (m-files)](#processing-of-data-from-hydrodynamic-codes-m-files)
  - [Data files (mat-files that can be loaded to workspace and used by Simulink templates)](#data-files-mat-files-that-can-be-loaded-to-workspace-and-used-by-simulink-templates)
  - [Hydrodynamics (m-files)](#hydrodynamics-m-files)
- [Frequency-domain identification (FDI) of radiation models (m-files)](#frequency-domain-identification-fdi-of-radiation-models-m-files)
  - [Demos](#demos)
  - [Utils](#utils)
---

## Simulink demos

```matlab
demoAUVdepthHeadingControl.slx          % simultaneously heading and depth control of the Remus 100 AUV
demoCS2passiveObserverDP.slx            %  passive observer with wave filtering and nonlinear PID control for the CyberShip2 model ship
demoDPThrusterModels.slx                % supply vessel with azimuth thrusters
demoDSRVdepthControl.slx                % depth control of DSRV
demoKalmanWavefilterAutop.slx           %  Kalman-filter based wave filter and heading autopilot for the mariner class cargo ship
demoNavalVesselMano.slx                 % zigzag test for the naval ship Mano  
demoNPSAUV.slx                          % NPS AUV heading control system
demoOtterUSVHeadingControl.slx	        % Otter USV heading control system
demoOtterUSVPathFollowingCourseControl  % Otter USV LOS path-following control using a course autopilot
demoOtterUSVPathFollowingHeadingControl % Otter USV ILOS and ALOS path-following control using a heading autopilot
demoPanamaxContainerShip.slx            % Panama container ship simulator
demoPassiveWavefilterAutopilot1.slx     % passive wave filter and heading autopilot design using compass measurements only
demoPassiveWavefilterAutopilot2.slx     % passive wave filter and heading autopilot design using compass and yaw rate measurements
demoS175WindCurrentAutopilot.slx        % S175 heading autopilot with wind and current loads
demoSemisubDPsystem.slx                 % semisubmersible DP system
demoWaveElevation.slx                   % computation of wave elevation from wave spectra
demoWaypointGuidance.slx                % waypoint guidance system
```

## GNC demos (m-files)
```matlab
GNCdemo        % MAIN PROGRAM
1) KinDemo     % Euler angle and quaternion kinematics
2) ManDemo     % Maneuvering trials
3) StabDemo    % Straight-line, directional and positional motion stability
4) WaveDemo    % Wave spectra demonstrations
```

## Examples (m-files)
```matlab
ExEKF           % for-loop implementation (predictor-corrector representation)  of a discrete-time extended Kalman filter (EKF) 
ExFeedback      % for-loop implementation for numerical integration of a 1st-order system under feedback and feedforward control
ExFFT           % estimation of the wave encounter frequency from time series using the fast-Fourier transform (FFT)
ExHybrid        % computation of a hybrid continuous path parametrized by waypoints
ExINS_AHRS      % Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and AHRS attitude measurements 
ExINS_Euler     % Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and compass measurements
ExINS_MEKF      % unit quaternion error-state (indirect) Kalman filter for INS aided by position and magnetic field measurements 
ExKF            % for-loop implementation (predictor-corrector representation) of a discrete-time linear Kalman filter (KF) 
ExKT            % computation of the Nomoto gain K and time constant T from a step response using nonlinear least-squares
ExLinspec       % linear approximations to the PM, JONSWAP, and Torsethaugen spectra using nonlinear least-squares
ExLQFinHor      % LQ finite time-horizon tracking controller for a mass-damper-spring system
ExLQR           % computes the LQR gains for a mass-damper system
ExLQtrack       % computes the LQ optimal tracking gains for a mass-damper system
ExMDS           % plots the step response of a 2nd-order mass-damper system
ExMSI           % plots the ISO 2631-1 (1997) and O'Hanlon and McCauley (1974) Motion Sickness Incidence (MSI) curves
ExNomoto        % Bode plots of ships parametrized by Nomoto’s 1st- and 2nd-order models
ExObsCtr        % observability and controllability matrices of a supply vessel
ExOtter         % simulates an Otter USV equipped with two propellers
ExPassiveObs    % plots the loop transfer function of the passive observer used for heading control
ExPathGen       % path generation using cubic polynomials 
ExPlotRAO       % script for plotting motion and force RAOs
ExPullout       % performs a pullout maneuver for two different ships
ExQuadProg      % quadratic programming applied to waypoint trajectory generation
ExQuest         % 6-DOF position/attitude vector from camera measurements using the QUEST algorithm
ExRefMod        % 2nd-order reference model with nonlinear damping and velocity saturation
ExResonance     % computes the closed-form responses in heave, roll, and pitch for a marine craft exposed to regular waves
ExRRD1          % roll and sway-yaw transfer functions for the Son and Nomoto container ship
ExRRD2          % rudder-roll damping (RRD) system for the Son and Nomoto container ship
ExRRD3          % inverse response in roll for the Son and Nomoto container ship due to a right-half-plane zero  (non-minimum phase)  
ExSMC           % integral sliding mode control (SMC) design for heading control
ExSpline        % path generation using cubic Hermite spline interpolation 
ExSTA           % adaptive-gain super twisting algorithm (STA) for heading control
ExTurnCircle    % generates the turning circle for two different ships
ExWageningen    % computes thrust and torque curves for a propeller using the Wageningen B-series data
ExWindForce     % plots the wind coefficients by Isherwoods (1972) 
ExZigZag        % generates zigzag maneuvers for two different ships and the Remus 100 AUV
```

## Marine craft simulator (m-files)

### Marine craft models

```matlab
clarke83        % linear maneuvering model parametrized using (L,B,T) found from linear regression of model tests (Clarke et al. 1983)
container       % nonlinear maneuvering model of a high-speed container ship, L = 175 m, including roll (Son and Nomoto 1982)
DSRV            % deep submergence rescue vehicle (DSRV), L = 5.0 m (Healey 1992)
frigate         % nonlinear autopilot model for a frigate, L = 100 m
Lcontainer      % linearized model of a high-speed container ship, L = 175 m, including the roll mode (Son and Nomoto 1982)
mariner         % nonlinear maneuvering model for the Mariner class vessel, L = 160 m 
navalvessel     % nonlinear maneuvering model of a multipurpose naval vessel, L = 51.5 m
npsauv          % Naval Postgraduate School autonomous underwater vehicle (AUV), L = 5.3 m 
otter           % OTTER small autonomous USV, L = 2.0 m
remus100        % REMUS 100 autonomous underwater vehicle (AUV), L = 1.9 m
rig             % semisubmersible linear mass-damper-spring model, L = 84.6 m
ROVzefakkel     % nonlinear autopilot model  of a boat, L = 45 m
supply          % linear DP model of a supply vessel, L = 76.2 m
tanker          % nonlinear course unstable maneuvering model of a tanker, L = 304 m 
```

### Time-series simulation

```matlab
SIMclarke83     % simulate clarke83.m under PD control
SIMmariner      % simulate mariner.m under PD control
SIMotter        % simulate otter.m under feedback control
SIMcontainer    % simulate container.m and Lcontainer.m under PD control
SIMnavalvessel  % simulate navalvessel.m under PD control
SIMremus100     % simulate remus100.m using autopilots for depth and heading control
SIMremus100ALOS % simulate remus100.m in a 3-D path-following scenario using adaptive LOS (ALOS)
SIMrig          % simulate the 6-DOF semisubmersible model under PID control
```

### Utils (m-files)

## Modelling

```matlab
addedMassSurge     % hydrodynamic added mass in surge, A11, approximated by the formula of Söding (1982)
coeffLiftDrag      % hydrodynamic lift and drag coefficients as a function of angle of attack of a submerged "wing profile"
forceLiftDrag      % hydrodynamic lift and drag forces as a function of angle of attack of a submerged "wing profile" 
forceSurgeDamping  % linear and quadratic damping forces in surge
crossFlowDrag 	   % crossflow drag computed from strip theory integrals
Dmtrx              % 6x6 linear damping matrix for marine craft (submerged and floating)
Gmtrx              % 6x6 system spring stiffness matrix G
gvect              % 6x1 vector of restoring forces, Euler angles as input
gRvect             % 6x1 vector of restoring forces, Euler angle, or unit quaternion rotation matrix as input
imlay61            % 6x6 hydrodynamic added mass and Coriolis-centripetal matrices MA and CA for a prolate spheroid
m2c                % 6x6 Coriolis-centripetal matrix C(nu) from system inertia matrix M
rbody              % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a general body
spheroid           % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a prolate spheroid 
wageningen         % thrust and torque coefficients of the Wageningen B-series propellers 
```

## Kinematics

```matlab
ecef2llh           % longitude, latitude, and height from ECEEF positions x, y, and z
euler2q            % unit quaternion from Euler angles
eulerang           % computes the Euler angle transformation matrices J, Rzyx and Tzyx
flat2llh           % longitude, latitude, and height from flat-earth positions x, y, and z
llh2ecef           % ECEEF positions x, y, and z from longitude, latitude, and height
llh2flat           % flat-earth positions x, y, and z from longitude, latitude, and height
R2euler            % Euler angles from rotation matrix elements
Rll                % Euler angle rotation matrix Rll for longitude and latitude
Rquat              % unit quaternion rotation matrix R in SO(3)
Rzyx               % Euler angle rotation matrix R in SO(3)
Tquat              % unit quaternion transformation matrix T, representing the attitude dynamics
Tzyx               % Euler angle transformation matrix T, representing the attitude dynamics
q2euler            % Euler angles from a unit quaternion
quatern            % unit quaternion transformation matrix J
quatprod           % quaternion product
quest              % quaternion rotation matrix R(q) and unit quaternion q between two vectors W = R(q) V
quest6dof          % 6-DOF vector eta = [x,y,z,phi,theta,psi] from three marker positions using the QUEST algorithm
```

## Environment

```matlab
blendermann94      % computes the wind forces and wind coefficients using Blendermann (1994)
encounter          % encounter frequency as a function of wave peak frequency, vessel speed, and wave direction
hs2vw              % converts significant wave height into an equivalent wind speed
isherwood72        % computes the wind forces and coefficients based on Isherwood (1972) 
rand_phases        % generates a uniformly distributed vector of random phases in the interval [-pi pi]
vw2hs              % converts average wind speed to significant wave height
waveresponse345    % steady-state heave, roll, and pitch responses for a ship in regular waves 
wavespec           % function used to evaluate different types of wave spectra
ww2we              % function used to transform a vector of wave frequencies to encounter frequencies
```

## Ship maneuvers

```matlab
pullout            % ship pullout maneuver
turncircle         % ship turning circle
zigzag             % zigzag maneuver for 3-DOF models
zigzag6dof         % zigzag maneuver for 6-DOF models
```

## Motion sickness

```matlab
ISOmsi             % ISO 2631-3, 1997 motion sickness incidence
HMmsi              % O'Hanlon and McCauley (1974) motion sickness incidence
```

## Transformations

```matlab
conversion         % defines global conversion factors for GNC applications
rad2pipi           % obsolete, use ssa
ssa                % smallest signed angle maps an angle in rad to the interval [-pi pi) or [-180 180)
Smtrx              % 3x3 vector skew-symmetric matrix S
Hmtrx              % 6x6 system transformation matrix H
vex                % computes a = vex(S(a)) where S is a skew-symmetric matrix
```

## Numerical integration methods

```matlab
euler2             % integrates a system of ordinary differential equations using Euler’s 2nd-order method
rk4                % integrates a system of ordinary differential equations using Runge-Kutta’s 4th-order method
```

### GNC (m-files)

## Guidance

```matlab
ALOS3D              % ALOS guidance laws for heading and pitch control in 3-D
ALOSpsi             % ALOS guidance law for heading control in 2-D (see demoOtterUSVPathFollowingHeadingControl.slx)
crosstrack          % computes the path-tangential origin and cross-track error for a target
crosstrackWpt       % computes the cross-track error when the path is a straight line between two waypoints
crosstrackWpt3D     % computes the 3-D tracking errors (along-, cross- and vertical-track errors)
hybridPath          % generates coefficients for sub-paths between waypoints
LOSchi              % LOS guidance law for course control in 2-D (see demoOtterUSVPathFollowingCourseControl.slx)
ILOSpsi             % ILOS guidance law for heading control in 2-D (see demoOtterUSVPathFollowingHeadingControl.slx)
order3              % path generation using cubic polynomials (see demoWaypointGuidance.slx)
order5              % path generation using 5th-order polynomials (see demoWaypointGuidance.slx)
```

## Navigation

```matlab
acc2rollpitch        % static roll and pitch angles from the specific force
EKF_5states          % estimation of SOG, COG, and course rate from NED positions or latitude-longitude
gravity              % acceleration of gravity as a function of latitude using the WGS-84 ellipsoid parameters
ins_ahrs             % error-state Kalman filter for INS aided by position and AHRS measurements 
ins_euler            % error-state Kalman filter for INS aided by position and yaw angle measurements
ins_mekf             % error-state Kalman filter for INS aided by position and magnetic field measurements
ins_mekf_psi         % error-state Kalman filter for INS aided by position and yaw angle measurements
insSignal            % basic INS signal generator
```

## Control

```matlab
lqtracker             % computes the LQ tracker gain matrices for LTI systems
nomoto                % generates Bode plots for the 1st- and 2nd-order Nomoto models
PIDnonlinearMIMO      % nonlinear MIMO PID regulator for dynamic positioning (DP)
ucalloc               % unconstrained control allocation
```

### HYDRO

## Hydrodynamic templates (Simulink)

```matlab
DP_ForceRAO.slx    % Simulink template for a DP vessel where wave loads are computed using force RAOs
DP_MotionRAO.slx   % Simulink template for a DP vessel where wave loads are computed using motion RAOs
MAN_ForceRAO.slx   % Simulink template for the unified maneuvering model where wave loads are computed using force RAOs
```

## Processing of data from hydrodynamic codes (m-files)

```matlab
veres2vessel       % reads data from ShipX output files and stores the data as a mat-file containing the structure <vessel>
vessel2ss          % computes the fluid-memory transfer functions and stores the data as a mat-file containing the structure <vesselABC>
wamit2vessel       % reads data from WAMIT output files and stores the data as a mat-file containing the structure <vessel>
```

## Data files (mat-files that can be loaded to workspace and used by Simulink templates)

```matlab
fpso, fpsoABC       % WAMIT data for a FPSO
semisub, semisubABC % WAMIT data for a semisubmersible
tanker, tankerABC   % WAMIT data for a tanker
s175, s175ABC       % ShipX data for a supply vessel
supply, supplyABC   % ShipX data for the S175
```

After loading the data files to the workspace using the Matlab command load, the following data structures are available:

| vessel   | vessel.main | vesselABC |
| -------- | --------    | --------  |
|          main: [1×1 struct]      |      name: 'tanker' | Ar: {6×6 cell} |
|    velocities: 1×60              |         T: draft    | Br: {6×6 cell} |
|      headings: [1×36 double]     |         B: beam     | Cr: {6×6 cell} |
|           MRB: [6×6 double]      |       Lpp: length   | Dr: {6×6 cell} | 
|             C: [6×6×60 double]   |         m: mass     | MRB: [6×6 double] | 
|         freqs: [1×60 double]     | rho: density of water | MA: [6×6 double] | 
|             A: [6×6×60 double]   | k44: radius of gyration | G: [6×6 double] | 
|             B: [6×6×60 double]   | k55: radius of gyration | Minv: [6×6 double] | 
|     motionRAO: [1×1 struct]      | k66: radius of gyration | r_g: [x_g y_g z_g] | 
|      forceRAO: [1×1 struct]      | g: acceleration of gravity | 
|      driftfrc: [1×1 struct]      | nabla: volume displacement | 
|            Bv: [6×6×60 double]   | CB: center of buyoancy |
|                                  | GM_T: tranverse metacentric height |
|                                  | GM_L: longitudinal metacentric height |
|                                  | CG: center of gravity | |


## Hydrodynamics (m-files)

```matlab
DPperiods           % periods and natural frequencies of a marine craft in DP
Hoerner             % 2-D Hoerner cross-flow form coefficient as a function of B and T
loadcond            % plots the roll and pitch periods as a function of GM_T and GM_L
plotABC             % plots the hydrodynamic coefficients Aij, Bij, and Cij as a function of frequency 
plotBv              % plots viscous damping Bvii as a function of frequency 
plotTF              % plots the motion or force RAO transfer functions
plotWD              % plots the wave drift amplitudes
```

## Frequency-domain identification (FDI) of radiation models (m-files)

### Demos
```matlab
Demo_FDIRadMod_NA     % FDI using hydrodynamic data without infinite-frequency added mass
Demo_FDIRadMod_WA     % FDI using hydrodynamic data including infinite-frequency added mass
```

### Utils
```matlab
EditAB                 % Function preparing the data for identification, select frequency range and remove wild-points
FDIRadMod              % Identify the SISO transfer function corresponding to the coupling specified
fit_siso_fresp         % Fit a continuous SISO transfer function to the frequency response data
ident_retardation_FD   % Identification of a parametric radiation convolution model K(s) = P(s)/Q(s)
ident_retardation_FDna % Identification of a parametric model A(jw) = B(w)/(jw) - A(w) 
```

