function [J,J1,J2] = quatern(q)
% [J,J1,J2] = QUATERN(q) computes the quaternion transformation matrices
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 6 October 2001, T I. Fossen - eta as first element in q 

eta  = q(1); eps1 = q(2); eps2 = q(3); eps3 = q(4); 
 
J1 = Rquat(q);
 
J2 = 0.5*[...
   -eps1 -eps2 -eps3        
    eta  -eps3  eps2
    eps3  eta  -eps1
   -eps2  eps1  eta   ];
 
J = [ J1  zeros(3,3);
      zeros(4,3) J2 ];
