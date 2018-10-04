% ExLQtrack  Computes the LQ optimal tracking gains for a mass-damper system
% Author:    Thor I. Fossen
% Date:      25 June 2001
% Revisions: 

% Design matrices
Q = diag([1]);      % user editable tracking error weights
R = diag([1]);      % user editable input weights
%
% System matrices
A = [0 1; -1 -2];   % user editable state matrix
B = [0; 1];         % user editable input matrix
C = [1 0];          % user editable output matrix

% Compute optimal gain matrices
[G1,G2] = lqtracker(A,B,C,Q,R)
