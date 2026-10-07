function [x,y,z] = llh2ecef(l,mu,h)
% [x,y,z] = LLH2ECEF(l,mu,h) computes the  ECEF positions (x,y,z)
% from longitude l (rad), latitude mu (rad) and height h above the surface 
% of the WGS-84 elipsoid.
%
% Author:    Thor I. Fossen
% Date:      2001-06-14
% Revisions:
%   2003-01-27 : Defined the l and mu inputs in radians.
%   2019-04-30 : Added decimals to r_e and r_p.
%   2026-10-07 : Updated the WGS-84 semi-minor axis (E. Krizman)

r_e = 6378137.0;            % WGS-84 data
r_p = 6356752.314245;

e2 = 1 - (r_p/r_e)^2;
N = r_e / sqrt(1 - e2 * sin(mu)^2);
x = (N + h) * cos(mu) * cos(l);
y = (N + h) * cos(mu) * sin(l);
z = (N * (1 - e2) + h) * sin(mu);
