function [l,mu,h] = flat2llh(xn,yn,zn,l0,mu0,h_ref)
% [l,mu,h] = flat2llh(xn,yn,zn,l0,mu0,h_ref) computes longitude l (rad), 
% latitude mu (rad) and height h (m) for the NED coordinates (xn,yn,zn) using
% a flat Earth coordinate system defined by the WGS-84 ellipsoid. The 
% flat Earth coordinate origin is located  at (l0, mu0) with reference 
% height h_ref in meters above the surface of the elipsoid. Both h and h_ref 
% are positive upwards, while z is postive downwards (NED).
%
% Author:    Thor I. Fossen
% Date:      20 July 2018
% Revisions: 2023-02-04 updatet the formulas for latitude and longitude
%
% Reference: 
% J. Farrell (2008). Aided Navigation: GPS with High Rate Sensors, 
% McGraw-Hill Professional, ISBN 9780071493291

% WGS-84 parameters
a = 6378137;             % Semi-major axis (equitorial radius)
f = 1 / 298.257223563;   % Flattening 
e = sqrt( 2*f - f^2 );   % Earth eccentricity

Rn = a / sqrt( 1-e^2 * sin(mu0)^2 );
Rm = Rn * ( ( 1-e^2 ) / ( 1-e^2 * sin(mu0)^2 ) );

% Equation (2.59) by Farrell (2008)
dmu = xn / ( Rm + h_ref);                   % delta latitude dmu = mu - mu0
dl  = yn / ( ( Rn + h_ref ) * cos(mu0) );   % delta longitude dl = l - l0

l = ssa(l0 + dl);
mu = ssa(mu0 + dmu);
h = h_ref - zn;
