% ExQuest    6 DOF position/attitude vector from camera measurements using the QUEST algorithm
%            See quest6DOF.m and quest.m     
% Author:    Thor I. Fossen and K. P. Lindegaard
% Date:      5 August 2001
% Revisions: 
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

% marker positions
mb1 = [0.69  0  -0.15]';
mb2 = [0.42  0  -0.51]';
mb3 = [-0.38 0  -0.18]';

% camera position
rcamera = [2.11  -0.48  1.89]';

% generate camera measurment for ship position: rcg = [ 2.5 3.0 0.1]' and 
% attitude: phi = 10 deg, theta = 5 deg and psi = 134 deg.
disp('True ship position and attitude:')
rcg = [ 2.5 3.0 0.1]'                        % loacation of CG
phi = 10
theta = 5
psi = 134
q = euler2q(phi*pi/180,theta*pi/180,psi*pi/180);  % quaternions
R = Rquat(q);                                     % rotation matrix
disp('camera measurments:')
y1 = rcg - rcamera + R*mb1                        % camera measurments
y2 = rcg - rcamera + R*mb2
y3 = rcg - rcamera + R*mb3

% QUEST ALGORITHM
% compute position/attitude eta in 6DOF, quaternions q and R from the camera measurements y1,y2,y3
y  = [y1; y2; y3];
mb = [mb1; mb2; mb3];
[eta,q,R] = quest6dof(y,mb,rcamera);
[phi,theta,psi] = q2euler(q);

disp('6 DOF position computed from camera measurements using the QUEST algorithm:')
rcg_quest    = eta(1:3)'
phi_quest    = phi*180/pi
theta_quest  = theta*180/pi
psi_quest    = psi*180/pi