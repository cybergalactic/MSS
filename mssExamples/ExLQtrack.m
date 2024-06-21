% ExLQtrack is compatible with MATLAB and GNU Octave (www.octave.org).
% The script computes the LQ optimal tracking gains for a mass-damper
% system.

% Author:    Thor I. Fossen
% Date:      2001-06-25
% Revisions: 

% Design matrices
Q = diag(1);      % User editable tracking error weights
R = diag(1);      % User editable input weights

% System matrices
A = [0 1; -1 -2];   % User editable state matrix
B = [0; 1];         % User editable input matrix
C = [1 0];          % User editable output matrix

% Compute optimal gain matrices
[G1,G2] = lqtracker(A,B,C,Q,R)
