function [xdot,U] = DSRV(x,u)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [xdot, U] = DSRV(in), with in=[x,u] returns returns the speed U in m/s
% (optionally) and the time derivative of the state vector 
% x = [ w q x z theta ]' for a Deep Submergence Rescue Vehicle (DSRV) of
% length 5.0 m, where
%
%   w:       heave velocity                 (m/s)
%   q:       pitch velocity                 (rad/s)
%   x:       x-position                     (m)
%   z:       z-position, positive downwards (m)
%   theta:   pitch angle                    (rad)
%
% Input:
%   u:       delta_s (rad), where delta_s is the stern plane
% 
% Reference: 
%   A. J. Healey (1992). Marine Vehicle Dynamics Lecture Notes and Problem 
%   Sets, Naval Postgraduate School (NPS), Monterey, CA.
% 
% Author:    Thor I. Fossen
% Date:      2001-05-12
% Revisions: 
%   2002-03-01 : Changed sign of Zdelta to positive.
%   2003-03-24 : Changed wrong speed 13.5 knots to 8 knots.
%   2021-06-21 : Converted to two input arguments [x, u]
%   2024-04-20 : Added compability to GNU Octave.

% Check of input and state dimensions
if (length(x)  ~= 5),error('x-vector must have dimension 5!'); end
if (length(u) ~= 1), error('u-vector must have dimension 1!'); end

% Cruise speed (m/s)       
U0 = 4.11;                   % U0 = 4.11 m/s = 8 knots = 13.5 ft/s
W0 = 0;      

% Normalization variables
L = 5.0;
U = sqrt( U0^2 + (W0 + x(1))^2 );

% states and inputs (with dimension)
delta_s = u; 
w       = x(1);  
q       = x(2);  
theta   = x(5);  
 
% Parameters, hydrodynamic derivatives and main dimensions
delta_max = deg2rad(30);         % max stern plane angle (deg)

Iy  =  0.001925;
m   =  0.036391;

Mqdot  = -0.001573; Zqdot  = -0.000130;
Mwdot  = -0.000146; Zwdot  = -0.031545;
Mq     = -0.01131;  Zq     = -0.017455;
Mw     =  0.011175; Zw     = -0.043938;
Mtheta = -0.156276 / U^2; 
Mdelta = -0.012797; Zdelta = 0.027695;

% Mass matrix elements
m11 = m-Zwdot;
m12 = -Zqdot;
m22 = Iy-Mqdot;
m21 = -Mwdot;

detM = m11 * m22 - m12 * m21;

% Rudder saturation
delta_s = sat(delta_s, delta_max);

% Forces and moments
Z = Zq * q + Zw * w                  + Zdelta * delta_s;
M = Mq * q + Mw * w + Mtheta * theta + Mdelta * delta_s;

% State derivatives (with dimension)
xdot = [ (m22 * Z - m12 * M) / detM
        (-m21 * Z + m11 * M) / detM
          cos(theta) * U0 + sin(theta) * w
         -sin(theta) * U0 + cos(theta) * w
              q                               ];
end
  