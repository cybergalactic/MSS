% eLQR is compatible with MATLAB and GNU Octave (www.octave.org).
% The script computes the LQR gains for a mass-damper system.
% 
% Author:    Thor I. Fossen
% Date:      2001-06-25
% Revisions: 

% Design matrices
Q = diag(1);      % User editable tracking error weights (dim m x m)
R = diag(1);      % User editable input weights (dim r x r)

% System matrices
A = [0 1; -1 -2];   % User editable state matrix (dim n x n)
B = [0; 1];         % User editable input matrix (dim n x r)
C = [1 0];          % User editable output matrix (dim m x n)
%
% Compute the optimal feedback gain matrix G
[K,P,E] = lqr(A,B,C'*Q*C,R);
G = -K
