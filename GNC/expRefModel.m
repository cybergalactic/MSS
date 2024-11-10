function x = expRefModel(t, t_max, x_start, x_final, alpha)
% expRefModel is compatible with MATLAB and GNU Octave (www.octave.org). 
% Computes the transition from x_start to x_final over time t_max using the 
% exponential curve 
% 
%    transition_factor = (1 - exp(-alpha * (t / t_max).^2));
%
% with convergence rate alpha. A target transition percentage of 95% is reached
% at t = t_max for alpha = 3.0, which is the mininmal recommended value.
%
% INPUTS:
%   t       - Current time
%   t_max   - Maximum time for the transition
%   x_start - Starting value
%   x_final - Final value
%   alpha   - Convergence rate, alpha >= 3.0
%
% Example usage: Change the depth of an AUV from z_start 0 to z_final = 10 
% meters in t_final = 20 seconds, and initiate the dive after t_dive = 40 seconds.
%
%   alpha = 5.0; 
%   z_start = 0;
%   z_final = 10;
%   t_final = 20;
%   t_dive = 40;
%   % Desired depth z_d at time t
%   z_d = expRefModel((t-t_dive), t_max, z_final, z_start, alpha); 
%
% Author: Thor I. Fossen
% Date: 2024-11-09
% Revisions:

if alpha < 3
    disp('The convergence rate alpha should be larger or equal to 3.0')
end

    % Ensure time is within bounds
if t > t_max
    t = t_max;
elseif t < 0
    t = 0;
end

% Exponential transition factor
transition_factor = (1 - exp(-alpha * (t / t_max).^2));

% Compute the current value x based on transition factor
x = x_start + (x_final - x_start) * transition_factor;

end