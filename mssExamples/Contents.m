% MSS Example Files 
%
% Ref.: T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. John Wiley & Sons Ltd., 2nd edition
%
% Examples: 
% ExEKF        - Discrte-time EKF using different measurement/sampling-time frequencies
% ExKF         - Discrete-time KF using different measurement/sampling-time frequencies
% ExFeedback   - Numerical integration of 1st-order system with feedback and feedforward control laws
% ExINS_AHRS   - INS error-state Kalman filter with AHRS
% ExINS_Euler  - INS error-state Kalman filter using Euler angles
% ExINS_MEKF   - INS error-state MEKF using unit quaternions
% ExFFT        - FFT used to compute the encounter frequency
% ExKT         - Computation of Nomoto gain and time const. using nonlinear least-squares (NLS)
% ExLinspec    - Linear approx. to PM, JONSWAP and Torsethaugen spectra using NLS
% ExLQFinHor   - LQ finite time-horizon tracking of mass-damper-spring system
% ExLQR        - Computes the LQR gains for a mass-damper system
% ExLQtrack    - Computes the LQ optimal tracking gains for a mass-damper system
% ExMDS        - Plots the step response of a 2nd-order mass-damper system 
% ExMSI        - Plots the ISO 2631-1 and O'Hanlon and McCauley Motion Sickness Incidence curves
% ExNomoto     - Generates the Nomoto Bode plots for a cargo ship and tanker
% ExObsCtr     - Observability and Controllability of Ships
% ExOtter      - Simulates the Otter USV equipped with two propeller inputs.
% ExPassiveObs - Plots the loop transfer function of the passive observer
% ExPathGen    - Path generation using cubic polynominals
% ExPlotRAO    - Script for plotting motion and force RAOs
% ExPullout    - Produces a pullout maneuver for the mariner class vessel and a container ship
% ExQuadProg   - Quadratic programming for way-point trajectory generation
% ExQuest      - 6 DOF position/attitude vector from camera measurements using the QUEST algorithm
% ExRefMod     - 2nd-order reference model with nonlinear damping and velocity saturation
% ExResonance  - Plots the amplitude as a function of resonant frequencies for ships
% EXRRD1       - Roll and sway-yaw transfer functions for the Son and Nomoto container ship
% ExRRD2       - Rudder-Roll Damping (RRD) system for the Son and Nomoto container ship
% ExRRD3       - Inverse response in roll due to right-half-plane zero (non-minimum phase)
% ExSMC        - Integral sliding mode control applied to AUV control
% ExSpline     - Cubic Hermite and spline interpolation of waypoints 
% ExSTA        - Supertwisting adaptive sliding model control applied to AUV control
% ExTurnCircle - Generates the turning circle for the mariner class vessel and a container ship
% ExWageningen - Computes thrust and torque for the Wageningen B-series propellers
% ExWindForce  - Plotting the wind forces and moment using Isherwoods (1972) formulas
% ExZigZag     - Produces a zigzag maneuver for the mariner class vessel and a container ship
