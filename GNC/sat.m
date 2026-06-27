function x = sat(x, x_max)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org)
% x = sat(x, x_max) saturates the input x at the specified maximum absolute 
% value x_max. If the inputs are vectors, saturation is performed elementwise.
%
% Symmetric saturation:
%   -x_max <= x <= x_max
%
% For asymmetric saturation limits, use:
%   x = satlim(x, x_min, x_max)
%
% Examples:
%   tau = sat(tau, tau_max);          % -tau_max <= tau <= tau_max
%   n   = satlim(n, n_min, n_max);    % n_min <= n <= n_max
%
% Inputs:
%   x:     Scalar or vector of numerical values.
%   x_max: Maximum allowable absolute value.
%
% Outputs:
%   x:     Saturated values.
%
% Author: Thor I. Fossen
% Date: 2024-04-20
% Revisions:
%   2024-06-05 : Extended to accept vectors as inputs.
%   2026-06-27 : Added reference to satlim() for asymmetric saturation.

% Check input
if any(x_max <= 0)
    error('x_max must contain positive elements.');
end

% Symmetric saturation
x = min(max(x, -x_max), x_max);

end
