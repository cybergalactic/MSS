% ExLQR      Computes the LQR gains for a mass-damper system
% Author:    Thor I. Fossen
% Date:      25 June 2001
% Revisions: 

% Design matrices
Q = diag([1]);      % user editable tracking error weights (dim m x m)
R = diag([1]);      % user editable input weights (dim r x r)
%
% System matrices
A = [0 1; -1 -2];   % user editable state matrix (dim n x n)
B = [0; 1];         % user editable input matrix (dim n x r)
C = [1 0];          % user editable output matrix (dim m x n)
%
% Compute the optimal feedback gain matrix G
[K,P,E] = lqr(A,B,C'*Q*C,R);
G = -K