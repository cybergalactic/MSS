function Ainv = invQR(A)
% invQR(A)is compatibel with MATLAB and GNU Octave (www.octave.org).
% The function Ainv = invQR(A) computes the matrix inverse of a strictly 
% positive matrix (A > 0) or positive definite matrix (A = A' > 0) using QR 
% matrix inversion. 
%
% Inputs:
%   A     : Strictly positive (A > 0) or positive definite matrix (A = A' > 0)
%                      
% Outputs:
%   Ainv  : Inverse of the matrix A
%
% Method:
%   The function utilizes the QR decomposition method to compute the inverse.
%   QR decomposition is preferred in this case due to its numerical stability 
%   over other methods, particularly for matrices that may not be perfectly 
%   conditioned. The steps are as follows:
%     1. Perform QR decomposition on the matrix A.
%     2. Use the decomposition to compute the inverse of A.
%   The QR decomposition factorizes matrix A into an orthogonal matrix Q and 
%   an upper triangular matrix R such that A = Q*R. Using the permutation 
%   matrix P, we can represent the inverse of A in a numerically stable way.
%
% Author: Thor I. Fossen
% Date: 2024-07-31
% Revisions:

% Check if A is near-singular
if rank(A) < size(A, 1)
    error('The matrix is near-singular and not invertible.');
end

% QR decomposition
[Q, R, P] = qr(A);

% Compute the inverse of A
Ainv = P * (R \ (Q'));

end