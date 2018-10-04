function [l,mu,h] = flat2llh(x,y,z,l0,mu0,h_ref)
% [l,mu,h] = FLAT2LLH(x,y,z,l0,mu0,h_ref) computes longitude l (rad), 
% latitude mu (rad) and height h (m) for the NED coordinates (x,y,z) using
% a flat Earth coordinate system defined by the WGS-84 ellipsoid.  
% The flat Earth coordinate origin is located  at (l0, mu0) with reference 
% height h_ref in meters above the surface of the Earth. Both h and h_ref 
% are positive upwards, while z is postive downwards (NED).
%
% Author:    Thor I. Fossen
% Date:      20 July 2018
% Revisions: 

% WGS-84 parameters
a = 6378137;           % Semi-major axis (equitorial radius)
f = 1/298.257223560;   % Flattening 
e = sqrt(2*f - f^2);   % Earth eccentricity

Rn = a / sqrt( 1-e^2 * sin(mu0)^2 );
Rm = Rn * ( (1-e^2) / (1-e^2 * sin(mu0)^2) );

dmu = x * atan2(1,Rm);
dl  = y * atan2(1,Rn*cos(mu0));

l = rad2pipi(l0 + dl);
mu = rad2pipi(mu0 + dmu);
h = h_ref - z;
