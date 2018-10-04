function q = euler2q(phi,theta,psi)
% q = EULER2Q(phi,theta,psi) computes the unit quaternions q = [eta eps1 eps2 eps3]
% from Euler angles phi, theta and psi
%
% Author:   Thor I. Fossen
% Date:     8th June 2000
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


R   = Rzyx(phi,theta,psi);
R(4,4) = trace(R);
[Rmax,i] = max( [R(1,1) R(2,2) R(3,3) R(4,4)] );
p_i= sqrt(1+2*R(i,i)-R(4,4));
if i==1,
   p1 = p_i;
   p2 = (R(2,1)+R(1,2))/p_i;
   p3 = (R(1,3)+R(3,1))/p_i;
   p4 = (R(3,2)-R(2,3))/p_i;
elseif i==2,
   p1 = (R(2,1)+R(1,2))/p_i;
   p2 = p_i;
   p3 = (R(3,2)+R(2,3))/p_i;
   p4 = (R(1,3)-R(3,1))/p_i;
elseif i==3,
   p1 = (R(1,3)+R(3,1))/p_i;
   p2 = (R(3,2)+R(2,3))/p_i;   
   p3 = p_i;
   p4 = (R(2,1)-R(1,2))/p_i;   
else
   p1 = (R(3,2)-R(2,3))/p_i;
   p2 = (R(1,3)-R(3,1))/p_i;
   p3 = (R(2,1)-R(1,2))/p_i;   
   p4 = p_i;
end

q = 0.5*[p4 p1 p2 p3]';
q = q/(q'*q);
