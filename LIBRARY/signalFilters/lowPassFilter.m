function [xf_next, y] = lowPassFilter(xf, u, w_n, h)
% lowPassFilter is compatible with MATLAB and GNU Octave (www.octave.org).
% This function implements a first-order low-pass filter using exact
% discretization. The filter smooths the input signal u[k] and updates
% the filtered output y[k] based on the natural frequencies w_n and the
% sampling time h. This version supports both scalar and vector inputs.
% The filter transfer function is:
%
%    h_lowPass(s) = w_n / (s + w_n)
%
% Inputs:
%   xf   - The previous filtered value (scalar or vector).
%   u    - The current input value to be filtered (scalar or vector).
%   w_n  - The natural frequencies of the low-pass filter (vector or 
%          diagonal matrix). If scalar, all filter frequencies are equal. 
%   h    - The sampling time in seconds (scalar).
%
% Output:
%   xf_next - Propagated filter state at time k+1 (scalar or vector)
%   y       - Filtered output at time k (scalar or vector)
%
% Author: Thor I. Fossen
% Date: 2024-04-26
% Revisions:

% Ensure w_n is a vector if a diagonal matrix is provided
if ~isvector(w_n)
    w_n = diag(w_n);
end

%  Filtered output at time k 
y = xf;

% Compute the filter coefficient
phi = exp(-h * w_n);

% Update the filtered value at time k+1
xf_next = phi .* xf + (1 - phi) .* u;

end