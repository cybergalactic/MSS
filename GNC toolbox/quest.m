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
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

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