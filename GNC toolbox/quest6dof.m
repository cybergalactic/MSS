function [eta,q,R] = quest6DOF(y,mb,rcamera)
% [eta,q,R] = quest6DOF(y,mb,rcamera)  computes the 6 DOF vector eta = [x,y,z,phi,theta,psi]
% and unit quaternion vector q = [eta,eps1,eps2,eps3] from three 3x1 marker position measurement 
% vectors stacked in y = [y1_e, y2_e, y3_e] using the QUEST algorithm (see quest.m). 
% The 3x1 body-frame marker position vectors are stacked in mb = [m1_b, m2_b, m3_b].
% 
% Inputs: 
%    y  = [x1,y1,z1,x2,y2,z2,x3,y3,z3] vector of 3 camera measurements (earth-frame)
%    mb = [xm1,ym1,zm1,xm2,ym2,zm2,xm3,ym3,zm3] vector of 3 marker positions (body frame)
%    rcamera = [x_c,y_c,z_c] vector of camera position (earth-frame)
%
% Outputs:
%    eta - 6x1 Earth-fixed positions and attitude (Euler angles)
%    q   - 4x1 vector of unit quaternions (optionally)
%    R   - 3x3 rotation matrix from BODY to NED (optionally) 
%
% Author:    Karl-Petter Lindegaard and Thor I. Fossen
% Date:      5th August 2001
% Revisions:
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

% check input dimensions
[mx,my] = size(mb);
[nx,ny] = size(y);  
if my>mx,help quest6DOF; error('mb must be a column vector of 3x1 position vectors'); end
if ny>nx,help quest6DOF; error(' y must be a column vector of 3x1 position vectors'); end
  
% Map inputs to 3x1 position vectors
m1_b = mb(1:3);  m2_b = mb(4:6);  m3_b = mb(7:9);
y1_e = y(1:3);   y2_e = y(4:6);   y3_e = y(7:9);

% Defining three position difference vectors for quest.m
W=[y1_e-y3_e y2_e-y3_e y1_e-y2_e];
V=[m1_b-m3_b m2_b-m3_b m1_b-m2_b];
   
% Call quest to find the rotation matrix R and unit quaternion q
[R,q] = quest(W,V);

% compute 6 DOF position/attitude vector
eta(1:3) = rcamera + (y1_e+y2_e+y3_e - R*(m1_b+m2_b+m3_b))/3;  % mean position
[eta(4) eta(5) eta(6)] = q2euler(q);                           % Euler angles
