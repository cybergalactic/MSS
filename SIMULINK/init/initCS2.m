%initCS2  The CybershipII 3DOF linear DP model
%date: 30.01.2003
%author: Dag Abel Sveen
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

%Overall length
LOA = 1.255;

%Center of gravity wrt. Aft Point
XG = 0.67;
YG = 0.0;
ZG = 0.0;

%Origo wrt. to CG
xg = 0.0425;

%Mass matrix
m = 23.8;
Iz = 1.76;
Xudot = -2.0;
Yvdot = -10.0;
Yrdot = -0.0; %!!!
Nrdot = -1.0;

M = [m-Xudot    0           0;
    0           m-Yvdot     m*xg-Yrdot;
    0           m*xg-Yrdot  Iz-Nrdot];

%Damping matrix
Xu = -2.0;
Yv = -7.0;
Yr = -0.1;
Nr = -0.5;

D = [ -Xu     0    0;
        0    -Yv  -Yr;
        0    -Yr  -Nr];

%thruster configuration data
%----------------------------
%Bow thruster wrt. CG
T3Lx = 0.465;
T3Ly = 0;
%Stbd thruster wrt. CG
T2Lx = -0.54;
T2Ly = 0.075;
%Port thruster wrt. CG
T1Lx = -0.54;
T1Ly = -0.075;
%Strb rudder wrt. CG
R2Lx = -0.595;
R2Ly = T2Ly;
%Port rudder wrt. CG
R1Lx = R2Lx;
R1Ly = T1Ly;
%Thruster coefficients
k3Tp = 1.84e-4;
k3Tn = 1.88e-4;

k2Tp = 3.74e-3; k1Tp = k2Tp;
k2Tn = 5.05e-3; k1Tn = k2Tn;

k1Ldelta1 = 0.927;  k2Ldelta1 = k1Ldelta1;      
k1Ldelta2 = -0.557; k2Ldelta2 = k1Ldelta2;
k1Ddelta1 = 0.079;  k2Ddelta1 = k1Ddelta1;
k1Ddelta2 = 0.615;  k2Ddelta2 = k1Ddelta2;
%velocity corrections
k1Ln = 2.1e-2;  k2Ln = k1Ln;
k1Dn = 9.64e-3; k2Dn = k1Dn;

%force allocation matrix
T = [1      1       0       0       0;
     0      0       1       1       1;
    -T1Ly  -T2Ly    T3Lx    R1Lx    R2Lx];
