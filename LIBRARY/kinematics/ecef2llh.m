function [l,mu,h] = ecef2llh(x,y,z)
% [l,mu,h] = ecef2llh(x,y,z) computes the longitude l (rad), latitude mu (rad)
% and height h (m) above the surface of the WGS-84 elipsoid from the
% ECEF positions (x,y,z).
%
% Author:   Thor I. Fossen
% Date:     2001-06-07
% Revisions:
%   2002-09-01 : Replaced atan2(y/x) by atan2(y,x).
%   2002-09-02 : Added the height output h.
%   2003-01-27 : Defined angle outputs in radians.
%   2020-02-19 : Added decimals to the WGS-84 parameters.
%   2026-01-08 : Added robust tan iteration, polar guard, and iteration limit.
%   2026-10-07 : Recomputed height from the final latitude iterate (E. Krizman)

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

% Recompute height using the final latitude iterate
c2 = 1 / (1 + t0^2);
s2 = t0^2 / (1 + t0^2);
N  = r_e^2 / sqrt( r_e^2 * c2 + r_p^2 * s2 );
h  = p * sqrt(1 + t0^2) - N;

