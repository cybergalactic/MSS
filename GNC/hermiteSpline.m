function [w_path, x_path, y_path, dx_path, dy_path, pi_h, ...
    pp_x, pp_y, N_horizon] = hermiteSpline(wpt, Umax, h)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org).
% hermiteSpline computes paths and path tangents using Hermite spline 
% interpolation. It uses the Akima algorithm for 1-D interpolation to 
% produce piecewise polynomials with continuous 1st-order derivatives (C1). 
% Compared to the spline algorithm, the Akima algorithm produces fewer 
% undulations and is better suited to deal with quick changes between 
% flat regions. Compared to the pchip algorithm, the Akima algorithm is 
% not as aggressively flattened and is therefore still able to deal with 
% oscillatory data. The function calculates a Hermite spline path through a 
% set of waypoints, estimating the numer of intervals (moving horizon) needed
% to find the crosstrack error for a vehicle moving at speed Umax. It 
% provides the path variables, the x and y coordinates along the path, 
% their derivatives, and the angle of the path tangent relative to North.
%
% Inputs:
%   wpt    - A struct containing fields `pos.x` and `pos.y` that are vectors 
%            of x and y coordinates of waypoints.
%   Umax   - Maximum speed along the path (m/s).
%   h      - Sampling interval for calculating path intervals (s).
%
% Outputs:
%   w_path  - A vector of path parameter values.
%   x_path  - X coordinates evaluated along the path.
%   y_path  - Y coordinates evaluated along the path.
%   dx_path - Derivatives of x coordinates along the path.
%   dy_path - Derivatives of y coordinates along the path.
%   pi_h    - Heading angles along the path relative to North (radians).
%   pp_x    - Piecewise polynomial for x coordinates.
%   pp_y    - Piecewise polynomial for y coordinates.
%   N_horizon - Number of intervals corresponding to a horizon of 1 second,
%               calculated based on maximum speed and path interval length.
%
% Example:
%   wpt.pos.x = [0, 100, 200];
%   wpt.pos.y = [0, 50, 100];
%   [w_path, x_path, y_path, dx_path, dy_path, pi_h, pp_x, pp_y, N_horizon] = 
%       hermiteSpline(wpt, 5, 0.1);
%
% Author:    Thor I. Fossen
% Date:      2024-04-21
% Revisions: 
%   2024-06-22 Replaced pchip with makima.

% Calculate path length from waypoints
pathLength = 0;
for i = 2:length(wpt.pos.x)
    deltaLength = sqrt((wpt.pos.x(i) - wpt.pos.x(i-1))^2 + ...
        (wpt.pos.y(i) - wpt.pos.y(i-1))^2);
    pathLength = pathLength + deltaLength;
end

% Calculate the time to traverse the path at maximum speed
time = pathLength / Umax;

% Determine the number of intervals based on the time and sampling interval
N_interval = floor(time / h) + 1; 
deltaPath = pathLength / N_interval;
N_horizon = round(Umax / deltaPath);

% Parameterize path using the calculated intervals
w_path = linspace(0, N_interval, N_interval + 1);
wpt.idx = linspace(0, N_interval, length(wpt.pos.x));

% Interpolate waypoints using MAKIMA
pp_x = makima(wpt.idx, wpt.pos.x);
pp_y = makima(wpt.idx, wpt.pos.y);

% Calculate derivatives of the path using the custom derivative function
pp_dx = ppDerivative(pp_x);
pp_dy = ppDerivative(pp_y);

% Calculate the heading angles relative to North
dx_path = ppval(pp_dx, w_path);
dy_path = ppval(pp_dy, w_path);
pi_h = atan2(dy_path, dx_path); % Heading angles in radians

% Evaluate the x and y coordinates at each path point
x_path = ppval(pp_x, w_path);
y_path = ppval(pp_y, w_path);

end

%% Calculate the derivative of a piecewise polynomial using finite differences
%function dpp = ppDerivative(pp)
%    [breaks, coefs, ~, k, ~] = unmkpp(pp);
 %   dcoefs = coefs(:, 1:k-1) .* (k-1:-1:1);
 %   dpp = mkpp(breaks, dcoefs);
% end

%% Calculate the derivative of a piecewise polynomial with smoothing
function dpp = ppDerivative(pp)
    % Extract the pieces of the piecewise polynomial
    [breaks, coefs, ~, k, ~] = unmkpp(pp);
    
    % Calculate the new coefficients for the derivative
    dcoefs = coefs(:, 1:k-1) .* (k-1:-1:1);
    dpp = mkpp(breaks, dcoefs);
    
    % Evaluate the polynomial and its derivative at a dense set of points
    num_points = 1000;
    xq = linspace(breaks(1), breaks(end), num_points);
    yq = ppval(pp, xq);
    dyq = ppval(dpp, xq);
    
    % Apply a simple moving average to smooth the derivatives
    window_size = 100; % Adjust this parameter to control smoothness
    smoothed_dyq = movmean(dyq, window_size);
    
    % Fit a new piecewise polynomial to the smoothed derivative
    dpp = interp1(xq, smoothed_dyq, 'linear', 'pp');

end


