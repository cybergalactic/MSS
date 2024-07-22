function [LOSangle, idx_start, y_e] = ...
    crosstrackHermiteLOS(w_path, x_path, y_path, dx_path, dy_path, pi_h, ...
    x, y, h, Delta_h, pp_x, pp_y, idx_start, N_horizon, gamma_h)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org).
% crosstrackHermiteLOS computes the Line-Of-Sight (LOS) angle for path 
% followingin a cubic Hermite spline path through a 2xn table of waypoints. 
% This function determines the desired course angle (chi_ref) or yaw angle 
% (psi_ref) using LOS guidance law.
%
% When gamma_h is omitted as a function argument, the desired course angle,
% chi_ref, is returned. If gamma_h is provided, the desired yaw angle, 
% psi_ref,is computed based on adaptive LOS guidance.
%
% Usage:
%   [chi_d, idx_start, y_e] = crosstrackHermiteLOS(w_path, x_path, y_path, 
%      dx_path, dy_path, pi_h, x, y, h, Delta_h, pp_x, pp_y, idx_start, 
%      N_horizon)
%
%   [psi_d, idx_start, y_e] = crosstrackHermiteLOS(w_path, x_path, y_path, 
%      dx_path, dy_path, pi_h, x, y, h, Delta_h, pp_x, pp_y, idx_start, 
%      N_horizon, gamma_h)
%
% Inputs:
%   w_path, x_path, y_path - Waypoints and their corresponding Hermite 
%                            spline coordinates on the path.
%   dx_path, dy_path       - Derivatives of the x and y coordinates.
%   pi_h                   - Path-tangential angle with respect to 
%                            North (rad).
%   x, y                   - Current vehicle position coordinates (m).
%   h                      - Sampling time (s).
%   Delta_h                - Look-ahead distance (m).
%   pp_x, pp_y             - Piecewise polynomial structures for x and y 
%                            coordinates.
%   idx_start              - Starting index for the moving horizon.
%   N_horizon              - Number of intervals within the horizon 
%                            (typically 1s).
%   gamma_h (optional)     - Positive adaptive gain constant for yaw angle 
%                            adjustment.
%
% Outputs:
%   LOSangle               - Desired course angle chi_d or yaw angle 
%                            psi_d (rad).
%   idx_start              - Updated start index of the moving horizon.
%   y_e                    - Cross-track error expressed in the path-
%                            tangential frame (m).
%
% The function utilizes LOS and ALOS guidance laws:
%   chi_ref[k] = pi_h - atan(y_e[k] / Delta_h)
%   psi_ref[k] = pi_h - beta_hat[k] - atan(y_e[k] / Delta_h)
%   beta_hat[k+1] = beta_hat[k] + h * gamma_h * Delta_h * ...
%                   y_e[k] / sqrt(Delta_h^2 + y_e[k]^2)
%
% Examples of use include autonomous vehicle path following and dynamic 
% positioning systems in marine and aerial vehicles.
%
% References:
%   T. I Fossen (2023). An Adaptive Line-of-sight (ALOS) Guidance Law 
%   for Path Following of Aircraft and Marine Craft. IEEE Transactions on 
%   Control Systems Technology 31(6),2887-2894. 
%   https://doi.org/10.1109/TCST.2023.3259819
%
% Author:    Thor I. Fossen
% Date:      2024-04-21
% Revisions:
%   None

persistent beta_hat  % Estimate of the crab angle

% Initialize beta_hat if empty
if isempty(beta_hat)
    beta_hat = 0;  % Initial parameter estimate
end

% Handle omitted gamma_h for course control
if nargin == 14
    gamma_h = [];
end

% Ensure index does not exceed bounds
idx_end = min(idx_start + N_horizon, length(w_path));
w_horizon = idx_start:idx_end;

% Evaluate polynomials over the horizon
x_horizon = ppval(pp_x, w_horizon);
y_horizon = ppval(pp_y, w_horizon);

% Calculate distances from vehicle position to all points on the horizon
distances = sqrt((x_horizon - x).^2 + (y_horizon - y).^2);

% Identify the minimum distance and corresponding index
[min_distance, min_distance_idx] = min(distances);

% Update the starting index for the moving horizon
idx_start = idx_start + min_distance_idx - 1;

% Calculate the signed cross-track error
vector_to_point = [(x - x_path(idx_start)), (y - y_path(idx_start))];
crossProd = dx_path(idx_start) * vector_to_point(2) ...
    - dy_path(idx_start) * vector_to_point(1);
y_e = sign(crossProd) * min_distance;

if nargin == 14  % LOS guidance law for course control
    LOSangle = pi_h(idx_start) - atan(y_e / Delta_h);
else  % ALOS guidance law for heading control
    LOSangle = pi_h(idx_start) - beta_hat - atan(y_e / Delta_h);
    % Update crab angle parameter estimate
    beta_hat = beta_hat + h * gamma_h * Delta_h * ...
               y_e / sqrt(Delta_h^2 + y_e^2);
end

end
