function x = sat(x, x_max)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org)
%
% x = sat(x, x_max) saturates the input x at the specified maximum absolute 
% value x_max.
%
% Inputs:
%   x:     A scalar, vector, or matrix of numerical values.
%   x_max: A positive scalar specifying the maximum allowable absolute value.
%
% Outputs:
%   x - The saturated values of the original input x
%
% Author: Thor I. Fossen
% Date: 2024-04-20

% Check if x_max is positive
if x_max < 0
    error('x_max must be a non-negative number.');
end

if abs(x) > x_max
    x = sign(x) * x_max;
end

end

