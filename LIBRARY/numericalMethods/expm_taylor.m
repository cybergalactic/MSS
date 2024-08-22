function E = expm_taylor(A)
% expm_taylor is compatible with MATLAB and GNU Octave (www.octave.org).
% The function computes the matrix exponential using a Taylor series 
% expansion. It serves as a replacement for the built-in expm function, 
% particularly in environments where a custom implementation is required.
%
% Inputs:
%   A - A square matrix of size n x n, for which the matrix exponential is 
%       to be computed.
%
% Outputs:
%   E - The matrix exponential of A, computed using a Taylor series expansion.
%
% Example: Matrix exponential of A
%   A = [0 1; -1 0];
%   E = expm_taylor(A);
%
% Example: INS unit quaternion propagation:
%   q_ins = expm_taylor( Tquat(w_ins) * h ) * q_ins; % Exact discretization
%   q_ins = q_ins / sqrt(q_ins' * q_ins);           % Normalization
%
% Notes:
%   - The function computes the matrix exponential by summing the terms of
%     the Taylor series until the norm of the term being added is smaller
%     than a tolerance (1e-6 in this case).
%   - For very large matrices or matrices with large elements, this method 
%     may converge slowly compared to more advanced techniques used in the 
%     built-in expm.
%
% See also:
%   expm, factorial

n = size(A, 1);  % Assume A is square
E = eye(n);  % Start with the identity matrix
term = eye(n);  % First term in the series (A^0/0!)
k = 1;
while norm(term, 'inf') > 1e-6 % Continue until terms are small
    term = (A^k) / factorial(k); % Compute next term
    E = E + term;  % Add term to sum
    k = k + 1;
end

end