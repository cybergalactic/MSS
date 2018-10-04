function [xdot,U] = DSRV(in)
% [xdot, U] = DSRV(in), with in=[x,ui] returns returns the speed U in m/s (optionally) and the 
% time derivative of the state vector x = [ w q x z theta ]' for a 
% deep submergence rescue vehicle (DSRV) L = 5.0 m, where
%
% w     = heave velocity                 (m/s)
% q     = pitch velocity                 (rad/s)
% x     = x-position                     (m)
% z     = z-position, positive downwards (m)
% theta = pitch angle                    (rad)
%
% The inputs are:
% ui     = delta (rad), where delta is the stern plane
% U0     = nominal speed (optionally). Default value is U0 = 4.11 m/s = 8 knots.
%
% Reference : A.J. Healey (1992). Marine Vehicle Dynamics Lecture Notes and 
%             Problem Sets, Naval Postgraduate School (NPS), Monterey, CA.
% 
% Author:    Thor I. Fossen
% Date:      12-May-2001
% Revisions: 01-Mar-2002: changed sign of Zdelta to positive
%            24-Mar-2003: changed typo 13.5 knots to 8 knots
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

% Check of input and state dimensions
x  = in(1:5);
ui = in(6);


if (length(x)  ~= 5),error('x-vector must have dimension 5 !'); end
if (length(ui) ~= 1),error('u-vector must have dimension 1 !'); end

% Cruise speed (m/s)       
U0 = 4.11;                   % U0 = 4.11 m/s = 8 knots = 13.5 ft/s
W0 = 0;      

% Normalization variables
L = 5.0;
U = sqrt( U0^2 + (W0+x(1))^2 );

% states and inputs (with dimension)
delta = ui; 

w     = x(1);  
q     = x(2);  
theta = x(5);  
 
% Parameters, hydrodynamic derivatives and main dimensions
delta_max = 30;             % max stern plane angle (deg)

Iy  =  0.001925;
m   =  0.036391;

Mqdot  = -0.001573; Zqdot  = -0.000130;
Mwdot  = -0.000146; Zwdot  = -0.031545;
Mq     = -0.01131;  Zq     = -0.017455;
Mw     =  0.011175; Zw     = -0.043938;
Mtheta = -0.156276/U^2; 
Mdelta = -0.012797; Zdelta = 0.027695;

% Masses and moments of inertia
m11 = m-Zwdot;
m12 = -Zqdot;
m22 = Iy-Mqdot;
m21 = -Mwdot;

detM = (m11*m22-m12*m21);

% Rudder saturation
if abs(delta) >= delta_max*pi/180, 
   delta = sign(delta)*delta_max*pi/180;
end

% Forces and moments
Z = Zq*q + Zw*w                + Zdelta*delta;
M = Mq*q + Mw*w + Mtheta*theta + Mdelta*delta;

% State derivatives (with dimension)
xdot = [ (m22*Z - m12*M)/detM
        (-m21*Z + m11*M)/detM
          cos(theta)*U0 + sin(theta)*w
         -sin(theta)*U0 + cos(theta)*w
              q                        ];