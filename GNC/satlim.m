function x = satlim(x, x_min, x_max)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org)
% x = satlim(x, x_min, x_max) saturates the input x between the specified lower 
% and upper limits. If the inputs are vectors, saturation is performed elementwise.
%
% Asymmetric saturation:
%   x_min <= x <= x_max
%
% For symmetric saturation about zero, use:
%   x = sat(x, x_max)
%
% Examples:
%   n   = satlim(n, n_min, n_max);
%   u   = satlim(u, 0, 1);
%   tau = sat(tau, tau_max);
%
% Inputs:
%   x:     Scalar or vector of numerical values.
%   x_min: Lower saturation limit.
%   x_max: Upper saturation limit.
%
% Outputs:
%   x:     Saturated values.
%
% Author: Thor I. Fossen
% Date: 2026-06-27

% Check input
if any(x_min >= x_max)
    error(['Each element of x_min must be strictly smaller than the ' ...
        'corresponding element of x_max.']);
end

% Asymmetric saturation
x = min(max(x, x_min), x_max);

end
