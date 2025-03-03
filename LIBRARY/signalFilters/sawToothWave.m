function y = sawToothWave(t)
% sawToothWave is compatible with MATLAB and GNU Octave (www.octave.org).
% The function generates a sawtooth wave signal which can be used for testing.
%
% Inputs:
%   t - time vector
%
% Inputs:
%   y - sawtooth signal vector
%
% Author: Thor I. Fossen
% Date: 2025-03-03
% Revision:

width = 0.5;

% Normalize time to the period (2*pi)
t_normalized = mod(t, 2*pi);

% Define the sawtooth wave based on width
y = zeros(size(t));

% For t < width * 2*pi, linear ramp up
increasing = t_normalized < (width * 2*pi);
y(increasing) = (t_normalized(increasing) / (width * 2*pi)) * 2 - 1;

% For t >= width * 2*pi, linear ramp down
decreasing = ~increasing;
y(decreasing) = -((t_normalized(decreasing) ...
    - width * 2*pi) / ((1-width) * 2*pi)) * 2 + 1;

end