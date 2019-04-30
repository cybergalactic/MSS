function [x,y,z] = llh2ecef(l,mu,h)
% [x,y,z] = LLH2ECEF(l,mu,h) computes the  ECEF positions (x,y,z)
% from longitude l (rad), latitude mu (rad) and height h
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 27th January 2003, inputs l and mu are defined in rad
%            30 Apr 2019, added decimals to r_e and r_p

r_e = 6378137; % WGS-84 data
r_p = 6356752.3142;

e = 0.08181979099211;
N = r_e^2/sqrt( (r_e*cos(mu))^2 + (r_p*sin(mu))^2 );
x = (N+h)*cos(mu)*cos(l);
y = (N+h)*cos(mu)*sin(l);
z = (N*(r_p/r_e)^2 + h)*sin(mu);