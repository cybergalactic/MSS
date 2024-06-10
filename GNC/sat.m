function x = sat(x, x_max)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org)
% x = sat(x, x_max) saturates the input x at the specified maximum absolute 
% value x_max. If the inputs are vectors, saturation is performed
% elementwise.
%
% Inputs:
%   x:     A vector of numerical values.
%   x_max: A vector specifying the maximum allowable absolute values.
%
% Outputs:
%   x:     A vector of saturated values of the original x elements
%
% Author: Thor I. Fossen
% Date: 2024-04-20
% Revisions:
%   2024-06-05 : Extended to accept vectors as inputs

% Check if x_max is positive
if any(x_max <= 0)
    error('x_max must be a vector of positive elements.');
end

% Saturation
x = min(max(x, -x_max), x_max);

end

