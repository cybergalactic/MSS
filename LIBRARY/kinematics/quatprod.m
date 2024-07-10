function q = quatprod(q1,q2)
% q = quatprod(q1,q2) computes the quaternion product
%
% Author:   Thor I. Fossen
% Date:     10th September 2010
% Revisions: 

eta1  = q1(1); 
eps1  = q1(2:4);
eta2  = q2(1); 
eps2  = q2(2:4);

q = [eta1*eta2-eps1'*eps2
     eta2*eps1+eta1*eps2+cross(eps1,eps2)];