function [phi,theta,psi] = q2euler(q)
% [phi,theta,psi] = Q2EULER(q) computes the Euler angles from the unit 
% quaternions q = [eta eps1 eps2 eps3]
%
% Author:   Thor I. Fossen
% Date:      2001-06-14  
% Revisions: 2007-09-03  Test for singular solution theta = +-90 deg has been improved
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
 

R = Rquat(q);
if abs(R(3,1))>1.0, error('solution is singular for theta = +- 90 degrees'); end
 
phi   = atan2(R(3,2),R(3,3));
theta = -asin(R(3,1));
psi   = atan2(R(2,1),R(1,1));
