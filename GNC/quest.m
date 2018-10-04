function [R,q] = quest(W,V)
% [R,q] = quest(W,V) computes the quaternion rotation matrix R(q) in SO(3)
% and the corresponding unit quaternion q = [eta eps1 eps2 eps3] (optionally output)
% between to vectors V and W such that: W = R(q) V. 
%
% Ref.: Schuster, M. D., Oh, S. D., 1981, Three-axis attitude determination from vector
%       observations, Journal of Guidance, Dynamics and Control JGC-4(1):70-77.
%
% Author:    Karl Petter Lindegaard
% Date:      20th June 2001
% Revisions: 1st August 2001 (Thor I. Fossen) - minor changes of notation

N = length(W(1,:));   % Number of column vectors
a = 1/N;

B = zeros(3,3);
for i=1:N,
    B = B + a*(W(:,i)*V(:,i)');
end

S   = B+B';
rho = trace(B);

Z = zeros(3,1);
for i=1:N,
    Z = Z + a*cross(W(:,i),V(:,i));
end

K = [ S-rho*eye(3) Z
      Z'           rho ];

[E,D] = eig(K);          % Find the eigenvector for largest eigenvalue of K
[a,b] = max(diag(D));

qq      = E(:,b);         % Compute the unit quaternion
qq      = qq/norm(qq);
q       = [qq(4) 
          -qq(1:3)];

R = Rquat(q);            % quaternion rotation matrix 