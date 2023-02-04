function [xn,yn,zn] = llh2flat(l,mu,h,l0,mu0,h_ref)
% [xn,yn,zn] = llh2flat(l,mu,h,l0,mu0,h_ref) computes (xn,yn,zn) for a flat Earth
% coordinate system from longitude l (rad), latitude mu (rad) and height h (m) 
% above the surface of the WGS-84 elipsoid. The flat Earth coordinate 
% origin is located  at (l0, mu0) with reference height h_ref in meters 
% above the surface of the WGS-84 elipsoid. Both h and h_ref are positive 
% upwards, while z is postive downwards (NED).
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

dl = l - l0;
dmu = mu - mu0;

Rn = a / sqrt( 1-e^2 * sin(mu0)^2 );
Rm = Rn * ( (1-e^2) / (1-e^2 * sin(mu0)^2) );

% Equation (2.59) by Farrell (2008)
xn = dmu * ( Rm + h_ref );
yn = dl * ( ( Rn + h_ref ) * cos(mu0) );
zn = h_ref - h;

