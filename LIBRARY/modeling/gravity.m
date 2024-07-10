function g = gravity(mu)
% g = gravity(mu) computes the acceleration of gravity (m/s^2) as a function
% of lattitude mu (rad) using the WGS-84 ellipsoid parameters.
% 
% Input:  mu  lattitude (rad)
%
% Author:    Thor I. Fossen
% Date:      19th February 2020
% Revisions: 

g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...
    sqrt( 1 - 0.006694384442 * sin(mu)^2 );
