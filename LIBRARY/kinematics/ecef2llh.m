function [l,mu,h] = ecef2llh(x,y,z)
% [l,mu,h] = ecef2llh(x,y,z) computes the longitude l (rad), latitude mu (rad)
% and height h (m) above the surface of the WGS-84 elipsoid from the
% ECEF positions (x,y,z).
%
% Author:   Thor I. Fossen
% Date:     7th June 2001
% Revisions: 1st September 2002, atan2(y/x) replaced by atan2(y,x)
%            2nd September 2002, new output argument for height h was added
%            27th January, 2003, angle outputs are defined in rad
%            19 February 2020, added decimals to WGS parameters
%            8 January 2026, robust tan-iteration + polar guard + max-iter

% WGS-84 data
r_e = 6378137.0;                
r_p = 6356752.314245;
e   = sqrt( 1 - (r_p/r_e)^2 );

% Longitude (four-quadrant)
l = atan2(y,x);

% Distance to spin axis
p = sqrt(x^2 + y^2);

% Polar guard 
p_eps = 1e-8;  % meters (safe tiny threshold)
if p < p_eps
    mu = sign(z) * pi/2;
    h  = abs(z) - r_p;
    return
end

% Iterate on t = tan(mu) (avoids sin/cos of mu during iteration)
tol   = 1e-10;
epsv  = 1;
k     = 0;
kmax  = 20;

% Initial guess t0 = tan(mu0) using spherical/ellipsoidal correction
t0 = (z/p) / (1 - e^2);

while (epsv > tol) && (k < kmax)

    % Compute cos^2(mu) and sin^2(mu) from t0
    c2 = 1 / (1 + t0^2);        % cos^2(mu)
    s2 = t0^2 / (1 + t0^2);     % sin^2(mu)

    % Prime vertical radius of curvature N (uses c2,s2 only)
    N  = r_e^2 / sqrt( r_e^2 * c2 + r_p^2 * s2 );

    % Height using 1/cos(mu) = sqrt(1+t^2)
    h  = p * sqrt(1 + t0^2) - N;

    % Fixed-point update for t = tan(mu)
    t  = (z/p) / ( 1 - e^2 * N/(N + h) );

    epsv = abs(t - t0);
    t0   = t;
    k    = k + 1;
end

% Latitude output (principal value in (-pi/2, pi/2))
mu = atan(t0);

