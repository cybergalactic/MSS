function E = expm_squaresPade(A)
    % Constants
    normA = norm(A, 'inf');  % Calculate the infinity norm of A
    maxDegree = 6;  % Degree of the Pade approximant, known to be stable/efficient

    % Scaling: Find the smallest s such that ||A/2^s|| <= 1/2
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

function [N, D] = pade_approx(A, m)
    % Compute the Pade approximant of degree (m, m)
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
