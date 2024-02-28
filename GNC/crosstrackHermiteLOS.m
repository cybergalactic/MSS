function [y_e, pi_h, chi_d, omega_chi_d, closestPoint, closestTangent] = ...
    crosstrackHermiteLOS(wayPoints, vehiclePos, Delta_h, U)
% [y_e, pi_h, chi_d, omega_chi_d, closestPoint, closestTangent] = ...
%    crosstrackHermiteLOS(wayPoints, vehiclePos, Delta_h, U) 
% computes the desired course angle when the path is a cubic Hermite spline 
% going through 2xn table of waypoints. The desired course angle chi_d and 
% course rate d/dt chi_d = omega_chi_d used by course autopilot systems are 
% computed using the proportional LOS guidance law:
%
%  chi_d = pi_h - atan( Kp * y_e ),    Kp = 1/Delta_h  
%
%  omega_chi_d = -Kp * U * sin( chi - pi_h ) / ( 1 + (Kp * y_e)^2 )
%
% The cross-track error y_e expressed in the path-tangential frame is
% computed as the minimum distance from the vehicle to the path tangent. If
% the desired angular rate omega_chi_d is to large for the vehicle, the
% speed U must be reduced during the course changing maneuver.
%
% Inputs:    
%  wayPoints: nx2 matrix of waypoints [x1, y1; x2, y2; ..., xn, yn] (m)
%  vehiclePos: [x, y] coordinate of the vehicle, used to compute the y_e
%  Delta_h: Look-ahead distance (m)
%  U: vehicle speed used to compute the desired course angle rate (rad/s)
%
% Outputs:
%  y_e: cross-track error expressed in the path-tangential frame (m)
%  pi_h: path-tangential angle with respect to NED (rad)
%  chi_d: desired LOS course angle (rad)
%  omega_chi_d: desired LOS course rate (rad/s)
%  closestPoint: closest point on the current segment to the vehicle
%  closestTangent: tangent at the closest point
%
% Calls: hermiteSpline.m
%
% See also: crosstrackWpt.m, LOSchi.m
%
% Ref. T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
% Motion Control. 2nd. Edition, Wiley
%
% Author:    Thor I. Fossen
% Date:      2024-02-28
% Revisions: 

n = size(wayPoints, 1);

% Initialize minimum distance to a large value
minDistance = inf;
closestPoint = [0, 0];
closestTangent = [0, 0];

% Add one dummy waypoint with bearing defined by the last two waypoints
R_inf = max(max(wayPoints));
bearing = atan2(wayPoints(end,2)-wayPoints(end-1,2), ...
    wayPoints(end,1)-wayPoints(end-1,1));
wayPoints(end+1,:) = wayPoints(end,:) + R_inf * [cos(bearing) sin(bearing)];

% Loop through each segment of the Hermite spline
for i = 1:n

    % Define the function to find the closest point on the current
    % segment to the vehicle
    distanceFunc = @(t) norm(hermiteSpline(t, i, wayPoints) - vehiclePos);

    % Find the t value that minimizes the distance for the current segment
    options = optimset('TolX',1e-6,'Display','off');
    t_min = fminbnd(distanceFunc, 0, 1, options);

    % Calculate the minimum distance for the current segment
    [pointOnSpline, tangentAtPoint] = hermiteSpline(t_min, i, wayPoints);
    minDist = norm(pointOnSpline - vehiclePos);

    if minDist < minDistance

        minDistance = minDist;

        % Update the closest point and tangent
        closestPoint = pointOnSpline;
        closestTangent = tangentAtPoint;

        % Determine the side using cross product (z-component)
        vectorToVehicle = vehiclePos - closestPoint;
        crossProd = tangentAtPoint(1) * vectorToVehicle(2)...
            - tangentAtPoint(2) * vectorToVehicle(1);

        % Positive if vehicle is on the left, negative if on the right
        side = sign(crossProd);

    end
end

% Signed cross-track error
y_e = minDistance * side;

% Desired course angle: chi_d (LOS guidance law) 
pi_h = atan2(closestTangent(2), closestTangent(1));
chi_d = pi_h - atan(y_e/Delta_h);

% Desired course rate: omega_chi_d
Dy_e = -U * y_e / sqrt( Delta_h^2 + y_e^2 );                   % d/dt y_e
omega_chi_d = -(1/Delta_h) * Dy_e / ( 1 + (y_e/Delta_h)^2 ); 


end