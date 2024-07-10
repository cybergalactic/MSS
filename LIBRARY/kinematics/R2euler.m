function [phi,theta,psi] = R2euler(R)
% [phi,theta,psi] = R2EULER(R) computes the Euler angles from the  
% the rotataion matirx S
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 

phi   = atan2(R(3,2),R(3,3));
theta = -asin(R(3,1));
psi   = atan2(R(2,1),R(1,1));