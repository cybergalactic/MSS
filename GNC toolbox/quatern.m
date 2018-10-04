function [J,J1,J2] = quatern(q)
% [J,J1,J2] = QUATERN(q) computes the quaternion transformation matrices
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 6 October 2001, T I. Fossen - eta as first element in q 
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

eta  = q(1); eps1 = q(2); eps2 = q(3); eps3 = q(4); 
 
J1 = Rquat(q);
 
J2 = 0.5*[...
   -eps1 -eps2 -eps3        
    eta  -eps3  eps2
    eps3  eta  -eps1
   -eps2  eps1  eta   ];
 
J = [ J1  zeros(3,3);
      zeros(4,3) J2 ];
