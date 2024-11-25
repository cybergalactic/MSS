function x = refModelPolyExp(t, t_max, x_start, x_final, method, alpha)
% refModelPolyExp is compatible with MATLAB and GNU Octave (www.octave.org). 
% Computes the transition from x_start to x_final over time t_max using
% one of three methods: 'exponential', 'cubic', or 'quintic'.
%
% INPUTS:
%   t       - Current time
%   t_max   - Maximum time for the transition
%   x_start - Starting value
%   x_final - Final value
%   method  - Transition method: 'exponential', 'cubic', or 'quintic'
%   alpha   - (Optional) Convergence rate for exponential method, alpha >= 3.0.
%             If not specified, the exponential method uses alpha = 3.0
%
% OUTPUT:
%   x       - The transitioned value at time t
%
% Example usage:
%   z_start = 0; z_final = 10; t_final = 20; t_dive = 40;
%   
%   Desired depth z_d at time t
%   z_d = refModelPolyExp((t - t_dive), t_final, z_start, z_final, 'exponential', 5); 
%   z_d = refModelPolyExp((t - t_dive), t_final, z_start, z_final, 'exponential'); 
%   z_d = refModelPolyExp((t - t_dive), t_final, z_start, z_final, 'qubic'); 
%   z_d = refModelPolyExp((t - t_dive), t_final, z_start, z_final, 'quintic'); 
%
% Author: Thor I. Fossen
% Date: 2024-11-23
% Revisions: 

% Ensure time is within bounds
if t > t_max
    t = t_max;
elseif t < 0
    t = 0;
end

% Normalized time
tau = t / t_max;

% Select transition method
switch method
    case 'exponential'
        % alpha is optional for exponential method
        if nargin < 6
            alpha = 3; 
        end

        % Validate alpha
        if alpha < 3
            error(['The convergence rate alpha should be larger or equal' ...
                ' to 3.0 for the exponential method.']);
        end
        % Exponential transition factor
        transition_factor = 1 - exp(-alpha * tau.^2);
        
    case 'cubic'
        % Cubic Hermite spline
        transition_factor = 3 * tau.^2 - 2 * tau.^3;
        
    case 'quintic'
        % Quintic Polynomial
        transition_factor = 10 * tau.^3 - 15 * tau.^4 + 6 * tau.^5;
        
    otherwise
        error(['Invalid method. Choose ''exponential'', ''cubic'', or ' ...
            '''quintic''.']);
end

% Compute the current value x based on transition factor
x = x_start + (x_final - x_start) * transition_factor;

end