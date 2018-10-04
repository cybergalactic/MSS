function g = gvect(W,B,theta,phi,r_g,r_b)
% g = GVECT(W,B,theta,phi,r_g,r_b) computes the 6x1 vector of restoring 
% forces about an arbitrarily point CO for a submerged body. For floating 
% vessels, see Gmtrx.m
% 
% Inputs:  W, B: weight and buoyancy
%          phi,theta: roll and pitch angles
%          r_g = [x_g y_g z_g]: location of CG with respect to CO
%          r_b = [x_b y_b z_b]: location of CB with respect to CO
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 20 oct 2008 improved documentation
%            22 sep 2013 corrected sign error on last row (Mohammad Khani)
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


sth  = sin(theta); cth  = cos(theta);
sphi = sin(phi);   cphi = cos(phi);

g = [...
   (W-B)*sth
  -(W-B)*cth*sphi
  -(W-B)*cth*cphi
  -(r_g(2)*W-r_b(2)*B)*cth*cphi + (r_g(3)*W-r_b(3)*B)*cth*sphi
   (r_g(3)*W-r_b(3)*B)*sth      + (r_g(1)*W-r_b(1)*B)*cth*cphi
  -(r_g(1)*W-r_b(1)*B)*cth*sphi - (r_g(2)*W-r_b(2)*B)*sth      ];
