function [LOSangle, LOSrate, y_e, pi_h, closestPoint, closestTangent] = ...
    crosstrackHermiteLOS(wayPoints, vehiclePos,h,U_h,Delta_h,gamma_h)
% [LOSangle, LOSrate, y_e, pi_h, closestPoint, closestTangent] = ...
%      crosstrackHermiteLOS(wayPoints,vehiclePos,h,U_h,Delta_h,gamma_h)
% computes the desired LOS angle when the path is a cubic Hermite spline 
% going through 2xn table of waypoints. The desired course angle and course 
% rate (chi_d, omega_chi_d) or desired yaw angle and yaw rate (psi_d, r_d) 
% used by course and heading utopilot systems are computed using the LOS or
% adaptive LOS (ALOS) guidance laws, respectively. The cross-track error 
% y_e expressed in the path-tangential frame is computed as the minimum 
% distance from the vehicle to the path tangent. If the desired angular 
% rate is to large for the vehicle, the horizontal speed U_h must be 
% reduced during the course changing maneuvers.
%
% Inputs:  
%  vehiclePos: Vehicle North-East position where
%    vehiclePos = [x y]
%  wayPoints = [wpt.pos.x wpt.pos.y] where
%    wpt.pos.x = [x1, x2,..,xn]' waypoints x-coord. expressed in NED (m)
%    wpt.pos.y = [y1, y2,..,yn]' waypoints y-coord. expressed in NED (m)
%  (x, y): Vehicle x-y coordinates (m), used to compute y_e
%  h: sampling time (s)
%  U_h: Vehicle speed (m/s) used to compute the desired turning rate
%  Delta_h: Look-ahead distance (m)
%  gamma_h: (OPTIONALLY) Positive adaptive gain constant. 
%
%  1) gamma_h unspecified and (chid_d, omega_chi_d) is computed using
%    [chi_d, omega_chi_d, y_e, pi_h, closestPoint, closestTangentt] = ...
%      crosstrackHermiteLOS(x,y,wpt,h,U_h,Delta_h)
%
%    LOS guidance law for course autopilot control: 
%      chi_d = pi_h - atan( Kp * y_e/Delta_h ),    Kp = 1/Delta_h  
%      omega_chi_d = -Kp * U * sin( chi - pi_h ) / ( 1 + (Kp * y_e)^2 )
%
%  2) gamma_h > 0 and (psi_d, r_d) is computed using 
%    [psi_d, r_d, y_e, pi_h, closestPoint, closestTangentt] = ...
%      crosstrackHermiteLOS(x,y,wpt,h,U_h,Delta_h,gamma_h)
%
%    Adaptive LOS (ALOS) guidance law for heading control:
%      psi_d = pi_h - beta_hat - atan( Kp * y_e ),    Kp = 1/Delta_h
%      d/dt beta_hat = gamma_h * Delta_h * y_e / sqrt( Delta^2 + y_e^2 )
%
%      d/dt y_e = -U * y_e / sqrt( Delta^2 + y_e^2 )
%      r_d = -Kp * dy_e/dt / ( 1 + (Kp * y_e)^2 ) 
%
% Outputs:
%  LOSangle: desired course angle chi_d or desired yaw angle psi_d (rad)
%  LOSrate: desired course rate omega_chi_d or desired yaw rate r_d (rad/s)
%  y_e: cross-track error expressed in the path-tangential frame (m)
%  pi_h: path-tangential angle with respect to NED (rad)
%  closestPoint: closest point on the current segment to the vehicle
%  closestTangent: tangent at the closest point
%
% Calls: hermiteSpline.m
%
% See also: crosstrackWpt.m, LOSchi.m, ALOSchi.m
%
% Ref. T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
% Motion Control. 2nd. Edition, Wiley
%
% Author:    Thor I. Fossen
% Date:      2024-03-29
% Revisions: 

persistent beta_hat;        % estimate of the crab angle

if isempty(beta_hat)
    beta_hat = 0;           % initial parameter estimate
end

if nargin == 5
    gamma_h = [];
end

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

% Guidance laws
pi_h = atan2(closestTangent(2), closestTangent(1));

if nargin == 5 % LOS guidance law 
    
    % Desired course angle: chi_d  
    LOSangle = pi_h - atan( y_e/Delta_h ); 

    % Desired course rate: omega_chi_d
    Dy_e = -U_h * y_e / sqrt( Delta_h^2 + y_e^2 ); % d/dt y_e
    LOSrate = -(1/Delta_h) * Dy_e / ( 1 + (y_e/Delta_h)^2 );

else % ALOS guidance law
    
    % Desired heading angle: psi_d  
    LOSangle = pi_h - beta_hat - atan( y_e/Delta_h );

    % Desired yaw rate: r_d
    Dy_e = -U_h * y_e / sqrt( Delta_h^2 + y_e^2 ); % d/dt y_e
    LOSrate = -(1/Delta_h) * Dy_e / ( 1 + (y_e/Delta_h)^2 );

    % integral state: y_int[k+1]
    beta_hat = beta_hat + ...
        h * gamma_h * Delta_h * y_e / sqrt( Delta_h^2 + y_e^2 );

end

end