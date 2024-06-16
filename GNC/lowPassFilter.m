function xf = lowPassFilter(xf, x, w_n, h)
% lowPassFilter is compatible with MATLAB and GNU Octave (www.octave.org). 
% This function implements a first-order low-pass filter using exact 
% discretization. The filter smooths the input signal `x` and updates
% the filtered output `xf` based on the natural frequency `w_n` and the 
% sampling time `h`.
%
% Inputs:
%   xf   - The previous filtered value (scalar).
%   x    - The current input value to be filtered (scalar).
%   w_n  - The natural frequency of the low-pass filter (rad/s).
%   h    - The sampling time (seconds).
%
% Outputs:
%   xf   - The updated filtered value (scalar).
%
% Example:
%   
%   w_n = 2 * pi * 1; % Natural frequency (1 Hz)
%   h = 0.01;         % Sampling period (0.01 seconds)
%   xf = 0;           % Initial filtered value
%
%   % Simulate input signal and apply low-pass filter
%   t = 0:h:10;       % Time vector
%   x = sin(2 * pi * 0.5 * t);  % Input signal (0.5 Hz sine wave)
%   xf_values = zeros(size(x)); % Preallocate output array
%
%   for i = 1:length(t)
%       xf = lowPassFilter(xf, x(i), w_n, h);
%       xf_values(i) = xf;
%   end
%
%   % Plot results
%   figure;
%   plot(t, x, 'b', t, xf_values, 'r');
%   xlabel('Time (s)');
%   ylabel('Amplitude');
%   legend('Input Signal', 'Filtered Signal');
%   title('Low-Pass Filter');
%
% Author: Thor I. Fossen
% Date: 2024-04-26
% Revisions:
%   None

phi = exp(-h * w_n);

xf = phi * xf + (1 - phi) * x;

end