function [J,J1,J2] = eulerang(phi,theta,psi)
% [J,Rzyx,Tzyx] = eulerang(phi,theta,psi) computes the Euler angle
% transformation matrix
%
%  J = [ Rzyx     0
%           0  Tzyx ]
%
% where J1 = Rzyx and J2 = Tzyx, see Rzyx.m and Tzyx.m.
%
% Author:    Thor I. Fossen
% Date:      14 Jun 2001
% Revisions: 08 May 2021  - added calls to Rzyx and Tzyx 
 
J1 = Rzyx(phi,theta,psi);
J2 = Tzyx(phi,theta);
 
J = [ J1  zeros(3,3)
      zeros(3,3) J2 ];