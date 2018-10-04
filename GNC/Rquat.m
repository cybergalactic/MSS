function R = Rquat(q)
% R = Rquat(q) computes the quaternion rotation matrix R in SO(3)
% for q = [eta eps1 eps2 eps3]
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 6 October 2001, T I. Fossen - eta as first element in q  

tol = 1e-6;
if abs(norm(q)-1)>tol; error('norm(q) must be equal to 1'); end

eta = q(1);
eps = q(2:4);

S = Smtrx(eps);
R = eye(3) + 2*eta*S + 2*S^2;