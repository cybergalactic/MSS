# MSS Quick Reference

A quick reference guide for the MATLAB MSS toolbox:

```matlab
>> help MSS               % List mss commands
>> help mssExamples       % List book examples (T. I. Fossen. Handbook of Marine Craft Hydrodynamics and Motion Control, Wiley 2021)
>> mssSimulink            % Simulink Library
```

## Table of Contents
- [Simulink demos](#simulink-demos)
- [Examples (m-files)](#examples-m-files)
- [Marine craft simulator (m-files)](#marine-craft-simulator-m-files)
  - [Marine craft models](#marine-craft-models)
  - [Time-series simulation](#time-series-simulation)
- [Modelling (m-files)](#modelling-m-files)
  - [Model components](#model-components)
  - [Kinematics](#kinematics)
  - [Environment](#environment)
  - [Transformations](#transformations)
  - [Numerical integration methods](#numerical-integration-methods)
- [HYDRO](#hydro)
  - [Hydrodynamic templates (Simulink)](#hydrodynamic-templatess-simulink)
  - [Processing of data from hydrodynamic codes (m-files)](#processing-of-data-from-hydrodynamic-codes-m-files)
  - [Data files (mat-files that can be loaded to workspace and used by Simulink templates)](#data-files-mat-files-that-can-be-loaded-to-workspace-and-used-by-simulink-templates)
- [GNC (m-files)](#gnc)
  - [Guidance](#guidance)
  - [Navigation](#navigation)
  - [Control](#control)
- [Frequency-domain identification (FDI) of radiation models (m-files)](#frequency-domain-identification-fdi-of-radiation-models-m-files)
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
ExRRD2          % rudder-roll damping (RRD) system for the Son and Nomoto  container ship
ExRRD3          % inverse response in roll for the Son and Nomoto container ship  due to a right-half-plane zero  (non-minimum phase)  
ExSMC           % integral sliding mode control (SMC) design for heading control
ExSpline        % path generation using cubic Hermite spline interpolation 
ExSTA           % adaptive-gain super twisting algorithm (STA) for heading control
ExTurnCircle    % generates the turning circle for two different ships
ExWageningen    % computes thrust and torque curves for a propeller using the Wageningen B-series data
ExWindForce     % plots the wind coefficients by Isherwoods (1972) 
ExZigZag        % generates zigzag maneuvers for two different ships
```

## Marine craft simulator (m-files)

### Marine craft models

```matlab
clarke83        % linear maneuvering model parametrized using (L,B,T) found from linear regression of model tests (Clarke et al. 1983)
container       % nonlinear maneuvering model of a high-speed container ship, L = 175 m, including the roll mode (Son and Nomoto 1982)
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
SIMremus100     % simulate remus100.m under PD control
SIMrig          % simulate the 6-DOF semisubmersible model under PID control
```

### Modelling (m-files)

## Model components

```matlab
addedMassSurge     % hydrodynamic added mass in surge, A11, approximated by the formula of Söding (1982)
coeffLiftDrag      % hydrodynamic lift and drag coefficients as a function of angle of attack of a submerged "wing profile"
forceLiftDrag      % hydrodynamic lift and drag forces as a function of angle of attack of a submerged "wing profile" 
forceSurgeDamping  % linear and quadratic damping forces in surge
crossFlowDrag 	   % crossflow drag computed from strip theory integrals
Dmtrx              % 6x6 linear damping matrix for marine craft (submerged and floating)
Gmtrx              % 6x6 system spring stiffness matrix G
gvect              % 6x1 vector of restoring forces, Euler angles as input
gRvect      			 % 6x1 vector of restoring forces, Euler angle, or unit quaternion rotation matrix as input
imlay61		         % 6x6 hydrodynamic added mass and Coriolis-centripetal matrices MA and CA for a prolate spheroid
m2c			           % 6x6 Coriolis-centripetal matrix C(nu) from system inertia matrix M
rbody              % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a general body
spheroid           % 6x6 rigid-body system inertia and Coriolis-centripetal matrices MRB and CRB of a prolate spheroid 
wageningen         % thrust and torque coefficients of the Wageningen B-series propellers 
```

## Kinematics
```matlab
```

## Environment
```matlab
```

## Transformations
```matlab
```

## Numerical integration methods
```matlab
```

### HYDRO

## Hydrodynamic templates (Simulink)

### Processing of data from hydrodynamic codes (m-files)
```matlab
```

### Data files (mat-files that can be loaded to workspace and used by Simulink templates)
```matlab
```

### GNC (m-files)

## Guidance
```matlab
```

## Navigation
```matlab
```

## Control
```matlab
```

### Frequency-domain identification (FDI) of radiation models (m-files)
```matlab
```

- 

