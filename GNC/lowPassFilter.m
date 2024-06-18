function xf = lowPassFilter(xf, x, w_n, h)
% lowPassFilter is compatible with MATLAB and GNU Octave (www.octave.org).
% This function implements a first-order low-pass filter using exact
% discretization. The filter smooths the input signal `x` and updates
% the filtered output `xf` based on the natural frequencies `w_n` and the
% sampling time `h`. This version supports both scalar and vector inputs.
%
% Inputs:
%   xf   - The previous filtered value (scalar or vector).
%   x    - The current input value to be filtered (scalar or vector).
%   w_n  - The natural frequencies of the low-pass filter (vector or 
%          diagonal matrix). If scalar, all filter frequencies are equal. 
%   h    - The sampling time in seconds (scalar).
%
% Outputs:
%   xf   - The updated filtered value (scalar or vector).
%
% Example:
%
%   w_n = 2 * pi * [1; 0.5; 0.2]; % Natural frequencies (1 Hz, 0.5 Hz, 0.2 Hz)
%   h = 0.01;         % Sampling period (0.01 seconds)
%   xf = [0; 0; 0];   % Initial filtered value for a 3-element vector
%
%   % Simulate input signal and apply low-pass filter
%   t = 0:h:10;            % Time vector
%   x = [sin(2 * pi * 0.5 * t); cos(2 * pi * 0.5 * t); t/10]; % 3-element vector input signal
%   xf_values = zeros(size(x)); % Preallocate output array
%
%   for i = 1:length(t)
%       xf = lowPassFilter(xf, x(:,i), w_n, h);
%       xf_values(:,i) = xf;
%   end
%
%   % Plot results
%   figure;
%   subplot(3,1,1);
%   plot(t, x(1,:), 'b', t, xf_values(1,:), 'r');
%   xlabel('Time (s)'); ylabel('Amplitude');
%   legend('Input Signal 1', 'Filtered Signal 1');
%   title('Low-Pass Filter for First Element');
%
%   subplot(3,1,2);
%   plot(t, x(2,:), 'b', t, xf_values(2,:), 'r');
%   xlabel('Time (s)'); ylabel('Amplitude');
%   legend('Input Signal 2', 'Filtered Signal 2');
%   title('Low-Pass Filter for Second Element');
%
%   subplot(3,1,3);
%   plot(t, x(3,:), 'b', t, xf_values(3,:), 'r');
%   xlabel('Time (s)'); ylabel('Amplitude');
%   legend('Input Signal 3', 'Filtered Signal 3');
%   title('Low-Pass Filter for Third Element');
%
% References:
%   None
%
% Author: Thor I. Fossen
% Date: 2024-04-26
% Revisions:
%   None

% Compute the filter coefficient
if isvector(w_n)
    phi = exp(-h * w_n);
else
    phi = diag(exp(-h * diag(w_n)));
end

% Update the filtered value
xf = phi .* xf + (1 - phi) .* x;

end