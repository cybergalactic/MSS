function E = expm_squaresPade(A)
% expm_squaresPade is compatible with MATLAB and GNU Octave (www.octave.org).
% The function computes the matrix exponential using the "Scaling and Squaring 
% method with Pade Approximants", which is the recommended method by Moler and
% Loan (2003). It serves as a replacement for the built-in expm function,
% particularly in environments where a custom implementation is required.
%
% Inputs:
%   A - A square matrix of size n x n, for which the matrix exponential is
%       to be computed.
%
% Outputs:
%   E - The matrix exponential of A.
%
% Example: Matrix exponential of A
%   A = [0 1; -1 0];
%   E = expm_squaresPade(A);
%
% Example: INS unit quaternion propagation:
%   q_ins = expm_squaresPade( Tquat(w_ins) * h ) * q_ins; % Exact discretization
%   q_ins = q_ins / norm(q_ins);                          % Normalization
%
% Reference:
%   C. Moler and C. V. Loan (2003). Nineteen Dubious Ways to Compute the 
%   Exponential of a Matrix, Twenty-Five Years Later. SIAM REVIEW, Vol. 45, No. 1.
%
% Author: Thor I. Fossen
% Date: 2024-08-28
% Revisions:

% Constants
normA = norm(A, 'inf');  % Calculate the infinity norm of A
maxDegree = 6;  % Degree of the Pade approximant, known to be stable/efficient

% Scaling: Find the smallest s such that || A/2^s || <= 1/2
s = max(0, ceil(log2(normA)) + 1);
A = A / 2^s;

% Pade Approximant of order (maxDegree, maxDegree)
[N, D] = pade_approx(A, maxDegree);
E = N / D;  % Approximation of e^(A/2^s)

% Squaring: Undo the scaling by repeated squaring
for i = 1:s
    E = E * E;
end

end

%% Compute the Pade approximant of degree (m, m)
function [N, D] = pade_approx(A, m)
% N: numerator matrix
% D: denominator matrix
I = eye(size(A));
X = A;
c = 1;
N = I;
D = I;
for k = 1:m
    c = c * (m + 1 - k) / (k * (2 * m + 1 - k));
    N = N + c * X;
    D = D + c * (-1)^k * X;
    X = A * X;  % Compute the next power of A
end

end
