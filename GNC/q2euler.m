function [phi,theta,psi] = q2euler(q)
% [phi,theta,psi] = Q2EULER(q) computes the Euler angles from the unit 
% quaternions q = [eta eps1 eps2 eps3]
%
% Author:   Thor I. Fossen
% Date:      2001-06-14  
% Revisions: 2007-09-03  Test for singular solution theta = +-90 deg has been improved

R = Rquat(q);
if abs(R(3,1))>1.0, error('solution is singular for theta = +- 90 degrees'); end
 
phi   = atan2(R(3,2),R(3,3));
theta = -asin(R(3,1));
psi   = atan2(R(2,1),R(1,1));
