function G = Gmtrx(nabla,A_wp,GMT,GML,r_p)
% G = GMTRX(nabla,A_wp,GMT,GML,r_gp) computes the 6x6 system spring stiffness matrix G
% about an arbitrarily point P for a floating vessel (small roll and pitch angles).
% For submerged vessels, see gvect.m
% 
% Inputs:  nabla: deplasement
%          Awp: water plane area
%          GMT, GML: transverse/longitudinal metacentric heights
%          r_p = [x_p y_p z_p]': location of P with respect to CO
%
% Author:     Thor I. Fossen
% Date:       14th June 2001
% Revisions:  26th June 2002,  variable Awp was replaced with A_wp 
%                              one zero in G was removed
%             20th Oct 2008,   improved documentation, use r_p for
%                              arbitrarily point
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


rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

Zz     = -rho*g*A_wp;
Kphi   = -rho*g*nabla*GMT;
Mtheta = -rho*g*nabla*GML;
G_CO   = diag([0 0 -Zz -Kphi -Mtheta 0]);  % assumes that CO = CF
G = Hmtrx(r_p)' * G_CO * Hmtrx(r_p);