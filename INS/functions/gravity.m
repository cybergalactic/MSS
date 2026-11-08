function g = gravity(mu)
% g = gravity(mu) computes the acceleration of gravity (m/s^2) as a function
% of latitude mu (rad) using the WGS-84 ellipsoid parameters.
% 
% Input:  mu  latitude (rad)
%
% Reference:
%   National Geospatial-Intelligence Agency (2014). Department of Defense
%   World Geodetic System 1984: Its Definition and Relationships with Local
%   Geodetic Systems. NGA.STND.0036_1.0.0_WGS84, Tables 3.5 and 3.6.
%
% Author:    Thor I. Fossen
% Date:      19th February 2020
% Revisions:
%   2026-11-08 : Updated WGS-84 normal gravity coefficients to the NGA (2014) values.

g = 9.7803253359 * ( 1 + 0.001931852652458 * sin(mu)^2 ) /...
    sqrt( 1 - 0.006694379990141 * sin(mu)^2 );
