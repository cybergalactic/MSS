function [J,J1,J2] = eulerang(phi,theta,psi)
% [J,J1,J2] = EULERANG(phi,theta,psi) computes the Euler angle
% transformation matrices
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 

cphi = cos(phi);
sphi = sin(phi);
cth  = cos(theta);
sth  = sin(theta);
 
J1 = Rzyx(phi,theta,psi);
 
if cth==0, error('J2 is singular for theta = +-90 degrees'); end
 
J2 = [...
      1  sphi*sth/cth  cphi*sth/cth;
      0  cphi          -sphi;
      0  sphi/cth      cphi/cth ];
 
J = [ J1  zeros(3,3);
      zeros(3,3) J2 ];