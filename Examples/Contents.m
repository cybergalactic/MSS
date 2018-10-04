% MSS Example Files 
%
% Ref.: T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. John Wiley & Sons Ltd.
%
% Examples
%    ExLQFinHor   - LQ finite time-horizon tracking of mass-damper-spring system
%    ExLQR        - Computes the LQR gains for a mass-damper system
%    ExLQtrack    - Computes the LQ optimal tracking gains for a mass-damper system
%    ExMDS        - Plots the step response of a 2nd-order mass-damper system 
%    ExKT         - Computation of Nomoto gain and time const. using nonlinear least-squares (NLS)
%    ExLinspec    - Linear approx. to PM, JONSWAP and Torsethaugen spectra using NLS
%    ExMSI        - Plots the ISO 2631-1 and O'Hanlon and McCauley Motion Sickness Incidence curves
%    ExNomoto     - Generates the Nomoto Bode plots for a cargo ship and tanker
%    ExObsCtr     - Observability and Controllability of Ships
%    ExPassiveObs - Plots the loop transfer function of the passive observer
%    ExPathGen    - Path generation using cubic polynominals
%    ExPullout    - Produces a pullout maneuver for the mariner class vessel and a container ship
%    ExQuadProg   - Quadratic programming for way-point trajectory generation
%    ExQuest      - 6 DOF position/attitude vector from camera measurements using the QUEST algorithm
%    ExRefMod     - 2nd-order reference model with nonlinear damping and velocity saturation
%    ExSpline     - Cubic Hermite and spline interpolation of waypoints 
%    EXRRD1       - Roll and sway-yaw transfer functions for the Son and Nomoto container ship
%    ExRRD2       - Rudder-Roll Damping (RRD) system for the Son and Nomoto container ship
%    ExRRD3       - Inverse response in roll due to right-half-plane zero (non-minimum phase)
%    ExTurnCircle - Generates the turning circle for the mariner class vessel and a container ship
%    ExWindForce  - Plotting the wind forces and moment using Isherwoods (1972) formulas
%    ExZigZag     - Produces a zigzag maneuver for the mariner class vessel and a container ship
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org
