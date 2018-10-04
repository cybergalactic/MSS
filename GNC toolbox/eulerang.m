function [J,J1,J2] = eulerang(phi,theta,psi)
% [J,J1,J2] = EULERANG(phi,theta,psi) computes the Euler angle
% transformation matrices
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
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


cphi = cos(phi);
sphi = sin(phi);
cth  = cos(theta);
sth  = sin(theta);
 
J1 = Rzyx(phi,theta,psi);
 
if cth==0, error('J2 is singular for theta = +-90 degrees'); end
 
J2 = [...
      1  sphi*sth/cth  cphi*sth/cth;
      0  cphi          -sphi;
      0  sphi/cth      cphi/cth ];
 
J = [ J1  zeros(3,3);
      zeros(3,3) J2 ];