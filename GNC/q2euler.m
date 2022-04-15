function [phi,theta,psi] = q2euler(q)
% [phi,theta,psi] = q2euler(q) computes the Euler angles from the unit 
% quaternion q = [eta eps1 eps2 eps3]
%
% Author:    Thor I. Fossen
% Date:      2001-06-14  
% Revisions: 2022-04-15 Fixed test to avoid abs(R(3,1)) > 1 and norm(q) > 1

q = q / norm(q);            % normalize q, handle round-off errors 
R = Rquat(q);

phi = atan2(R(3,2),R(3,3));

if (abs( R(3,1 )) > 1)      % handle NaN due to round-off errors
    R(3,1) = sign(R(3,1));
end

theta = -asin(R(3,1));
psi = atan2(R(2,1),R(1,1));
