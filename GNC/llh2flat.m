function [x,y,z] = llh2flat(l,mu,h,l0,mu0,h_ref)
% [x,y,z] = LLH2FLAT(l,mu,h,l0,mu0,h_ref) computes (x,y,z) for a flat Earth
% coordinate system from longitude l (rad), latitude mu (rad) and height h (m) 
% using the WGS-84 ellipsoid.  The flat Earth coordinate origin is located 
% at (l0, mu0) with reference height h_ref in meters above the surface of 
% the Earth. Both h and h_ref are positive upwards, while z is postive
% downwards (NED).
%
% Author:    Thor I. Fossen
% Date:      20 July 2018
% Revisions: 

% WGS-84 parameters
a = 6378137;           % Semi-major axis (equitorial radius)
f = 1/298.257223560;   % Flattening 
e = sqrt(2*f - f^2);   % Earth eccentricity

dl = l - l0;
dmu = mu - mu0;

Rn = a / sqrt( 1-e^2 * sin(mu0)^2 );
Rm = Rn * ( (1-e^2) / (1-e^2 * sin(mu0)^2) );

x = dmu / atan2(1,Rm);
y = dl / atan2(1,Rn*cos(mu0));
z = h_ref - h;
