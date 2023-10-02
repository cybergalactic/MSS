# MSS Quick Reference

A quick reference guide for MATLAB syntax and operations.

## Table of Contents
- [Simulink demos](#simulink-demos)
- [Examples (m-files)](#examples-mfiles)

---

mssSimulink.slx				 Simulink Library

## Simulink demos

```matlab
demoAUVdepthHeadingControl 		          % simultaneously heading and depth control of the Remus 100 AUV
demoCS2passiveObserverDP.slx		        %  passive observer with wave filtering and nonlinear PID control for the CyberShip2 model ship
demoDPThrusterModels.slx			          % supply vessel with azimuth thrusters
demoDSRVdepthControl.slx			          % depth control of DSRV
demoKalmanWavefilterAutop.slx		        %  Kalman-filter based wave filter and heading autopilot for the mariner class cargo ship
demoNavalVesselMano.slx			            % zigzag test for the naval ship Mano  
demoNPSAUV.slx				                  % NPS AUV heading control system
demoOtterUSVHeadingControl.slx	        % Otter USV heading control system
demoOtterUSVPathFollowingCourseControl  % Otter USV LOS path-following control using a course autopilot
demoOtterUSVPathFollowingHeadingControl % Otter USV ILOS and ALOS path-following control using a heading autopilot
demoPanamaxContainerShip.slx		        % Panama container ship simulator
demoPassiveWavefilterAutopilot1.slx	    % passive wave filter and heading autopilot design using compass measurements only
demoPassiveWavefilterAutopilot2.slx	    % passive wave filter and heading autopilot design using compass and yaw rate measurements
demoS175WindCurrentAutopilot.slx		    % S175 heading autopilot with wind and current loads
demoSemisubDPsystem.slx			            % semisubmersible DP system
demoWaveElevation.slx		   	            % computation of wave elevation from wave spectra
demoWaypointGuidance.slx			          % waypoint guidance system
```

## Examples (m-files)
```matlab
ExEKF			        % for-loop implementation (predictor-corrector representation)  of a discrete-time extended Kalman filter (EKF) 
ExFeedback		    % for-loop implementation for numerical integration of a 1st-order system under feedback and feedforward control
ExFFT 			estimation of the wave encounter frequency from time-series using the fast-Fourier transform (FFT)
ExHybrid		computation of a hybrid continuous path parametrized by waypoints
ExINS_AHRS 		Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and AHRS attitude measurements 
ExINS_Euler 		Euler angle error-state (indirect) Kalman filter for INS aided by GNSS position and compass measurements 
ExINS_MEKF 		unit quaternion error-state (indirect) Kalman filter for INS aided by position and magnetic field measurements 
ExKF			for-loop implementation (predictor-corrector representation)  of a discrete-time linear Kalman filter (KF) 
ExKT			computation of the Nomoto gain K and time constant T from a step response using nonlinear least-squares
ExLinspec		linear approximations to the PM, JONSWAP and Torsethaugen spectra using nonlinear least-squares
ExLQFinHor		LQ finite time-horizon tracking controller for a mass-damper-spring system
ExLQR 			computes the LQR gains for a mass-damper system
ExLQtrack		computes the LQ optimal tracking gains for a mass-damper system
ExMDS			plots the step response of a 2nd-order mass-damper system
ExMSI			plots the ISO 2631-1 (1997) and O'Hanlon and McCauley (1974) Motion Sickness Incidence (MSI) curves
ExNomoto		Bode plots of ships parametrized by Nomoto’s 1st- and 2nd-order models
ExObsCtr		observability and controllability matrices of a supply vessel
ExOtter		simulates an Otter USV equipped with two propellers
ExPassiveObs 		plots the loop transfer function of the passive observer used for heading control
ExPathGen 		path generation using cubic polynomials 
ExPlotRAO 		script for plotting motion and force RAOs
ExPullout 		performs a pullout maneuver for two different ships
ExQuadProg		quadratic programming applied to waypoint trajectory generation
ExQuest 		6-DOF position/attitude vector from camera measurements using the QUEST algorithm
ExRefMod 		2nd-order reference model with nonlinear damping and velocity saturation
ExResonance 		computes the closed-form responses in heave, roll and pitch for a marine craft exposed to regular waves
ExRRD1 		roll and sway-yaw transfer functions for the Son and Nomoto container ship
ExRRD2 		rudder-roll damping (RRD) system for the Son and Nomoto  container ship
ExRRD3 		inverse response in roll for the Son and Nomoto container ship  due to a right-half-plane zero  (non-minimum phase)  
ExSMC 			integral sliding mode control (SMC) design for heading control
ExSpline 		path generation using cubic Hermite spline interpolation 
ExSTA 			adaptive-gain super twisting algorithm (STA) for heading control
ExTurnCircle 		generates the turning circle for two different ships
ExWageningen 	computes thrust and torque curves for a propeller using the Wageningen B-series data
ExWindForce		plots the wind coefficients by Isherwoods (1972) 
ExZigZag		generates zigzag maneuvers for two different ships

```

