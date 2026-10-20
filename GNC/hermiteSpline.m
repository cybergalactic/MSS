function [w_path, x_path, y_path, dx_path, dy_path, pi_h, ...
    pp_x, pp_y, N_forward,kappa] = hermiteSpline(wpt, Umax, h)
% hermiteSpline computes paths and path tangents using Hermite spline 
% interpolation with chord-length parameterization. Interior tangents bisect
% the incoming and outgoing segment directions, giving a shared nonzero
% derivative through corners, including axis-aligned right-angle turns.
% Endpoint tangents follow the first and last segments. Consecutive duplicate
% waypoints and exact direction reversals are not supported.
% The function calculates a path through the waypoints, estimating the
% number of intervals (moving horizon) needed
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
%   x_path  - X coordinates evaluated along the path (North).
%   y_path  - Y coordinates evaluated along the path (East).
%   dx_path - Derivatives of x coordinates along the path.
%   dy_path - Derivatives of y coordinates along the path.
%   pi_h    - Path tangent angles along the path relative to North (radians).
%   pp_x    - Piecewise polynomial for x coordinates.
%   pp_y    - Piecewise polynomial for y coordinates.
%   N_forward - Number of intervals corresponding to a horizon of 1 second,
%               calculated based on maximum speed and path interval length.
%   kappa   - Curvature along the path.
%
% Example:
%   wpt.pos.x = [0, 100, 200];
%   wpt.pos.y = [0, 50, 100];
%   [w_path, x_path, y_path, dx_path, dy_path, pi_h, pp_x, pp_y, N_forward] = 
%       hermiteSpline(wpt, 5, 0.1);
%
% Author:    Thor I. Fossen
% Date:      2024-04-21
% Revisions: 
%   2024-06-22 Replaced pchip with makima.
%   2026-10-18 Replaced coordinate-wise MAKIMA with shared Hermite tangents.

validateattributes(Umax, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(wpt.pos.x, {'numeric'}, {'real','finite','vector','numel',numel(wpt.pos.y)});
validateattributes(wpt.pos.y, {'numeric'}, {'real','finite','vector'});
points = [wpt.pos.x(:), wpt.pos.y(:)];
if size(points, 1) < 2
    error('hermiteSpline:Waypoints', 'At least two waypoints are required.');
end
segments = diff(points, 1, 1);
lengths = hypot(segments(:,1), segments(:,2));
if any(lengths == 0)
    error('hermiteSpline:DuplicateWaypoints', ...
        'Consecutive waypoints must be distinct.');
end
directions = bsxfun(@rdivide, segments, lengths);
tangents = [directions(1,:); ...
    directions(1:end-1,:) + directions(2:end,:); directions(end,:)];
tangentNorms = hypot(tangents(:,1), tangents(:,2));
if any(tangentNorms < 1e-12)
    error('hermiteSpline:DirectionReversal', ...
        'Direction reversals require additional waypoints to define a turn.');
end
tangents = bsxfun(@rdivide, tangents, tangentNorms);
pathLength = sum(lengths);

% Calculate the time to traverse the path at maximum speed
time = pathLength / Umax;

% Determine the number of intervals based on the time and sampling interval
N_interval = floor(time / h) + 1; 
deltaPath = pathLength / N_interval;
N_forward = round(Umax / deltaPath);

% Retain integer sample parameters for the moving-horizon guidance search.
% Polynomial breaks are spaced in proportion to waypoint chord lengths.
w_path = linspace(0, N_interval, N_interval + 1);
breaks = [0; cumsum(lengths)] * (N_interval / pathLength);
tangents = tangents * (pathLength / N_interval);

% Construct cubic Hermite polynomials with the same tangent on both sides
% of every waypoint. Coefficients use the local parameter w - breaks(i).
widths = diff(breaks);
slopes = bsxfun(@rdivide, segments, widths);
cubic = bsxfun(@rdivide, ...
    tangents(1:end-1,:) + tangents(2:end,:) - 2 * slopes, widths.^2);
quadratic = bsxfun(@rdivide, ...
    3 * slopes - 2 * tangents(1:end-1,:) - tangents(2:end,:), widths);
pp_x = mkpp(breaks', [cubic(:,1), quadratic(:,1), ...
    tangents(1:end-1,1), points(1:end-1,1)]);
pp_y = mkpp(breaks', [cubic(:,2), quadratic(:,2), ...
    tangents(1:end-1,2), points(1:end-1,2)]);

% Calculate derivatives of the path using the custom derivative function
pp_dx = ppDerivative(pp_x);
pp_dy = ppDerivative(pp_y);

% Calculate the path tangent relative to North
dx_path = ppval(pp_dx, w_path);
dy_path = ppval(pp_dy, w_path);
pi_h = atan2(dy_path, dx_path); % Path tangent angles in radians

% Second derivatives for curvature computation
pp_ddx = ppDerivative(pp_dx);
pp_ddy = ppDerivative(pp_dy);

ddx_path = ppval(pp_ddx, w_path);
ddy_path = ppval(pp_ddy, w_path);

% Curvature of parametric path
kappa = (dx_path .* ddy_path - dy_path .* ddx_path) ./ ...
        (dx_path.^2 + dy_path.^2).^(3/2);

% Evaluate the x and y coordinates at each path point
x_path = ppval(pp_x, w_path);
y_path = ppval(pp_y, w_path);

end

%% Differentiate the piecewise polynomial coefficients analytically
function dpp = ppDerivative(pp)
    [breaks, coefs, ~, k, ~] = unmkpp(pp);
    dcoefs = coefs(:,1:k-1) .* (k-1:-1:1);
    dpp = mkpp(breaks, dcoefs);
end



