function [LOSangle, y_e, pi_h, closestPoint, closestTangent] = ...
    crosstrackHermiteLOS(wayPoints, vehiclePos, h, undulation_factor,...
    Delta_h, gamma_h)
% crosstrackHermiteLOS: computes the desired LOS angle when the path is a 
% cubic Hermite spline going through 2xn table of waypoints. The desired 
% course angle, chi_ref, is returned if gamma_h is omitted as function argument: 
%
%  [chi_d, y_e, pi_h, closestPoint, closestTangent] = ...
%     crosstrackHermiteLOS(wayPoints,vehiclePos,h,undulation_factor,Delta_h)
%
% Alternatively, by specifying gamma_h > 0, the desired yaw angle, psi_ref, 
% is computed:
%
%  [psi_d, y_e, pi_h, closestPoint, closestTangent] = ...
%     crosstrackHermiteLOS(wayPoints, vehiclePos, h, undulation_factor,...
%     Delta_h, gamma_h)
%
% The desired course angle is computed using the LOS guidance las
%
%  chi_ref[k] = pi_h - atan( y_e[k] / Delta_h )
%
% while the desired heading angle is computed using the adaptive LOS (ALOS) 
% guidance law by Fossen (2023),
%
%  psi_ref[k] = pi_h - beta_hat[k] - atan( y_e[k] / Delta_h )
%  beta_hat[k+1] = beta_hat[k] + h * gamma_h * Delta_h * ...
%      y_e[k] / sqrt( Delta_h^2 + y_e[k]^2 )
%
% The cross-track error, y_e, expressed in the path-tangential frame and 
% it is computed as the minimum distance from the vehicle to the path 
% tangent. To handle steps in the course and heading commands due to 
% waypoint switching, the following observer can be applied:
%
%  [chi_d, omega_chi_d] = LOSobserver(chi_d, omega_chi_d, chi_ref, h, K_f, T_f)
%  [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f, T_f)
%
% Inputs:  
%
%  vehiclePos: Vehicle North-East position where
%    vehiclePos = [x y]
%  wayPoints = [wpt.pos.x wpt.pos.y] where
%    wpt.pos.x = [x1, x2,..,xn]' waypoints x-coord. expressed in NED (m)
%    wpt.pos.y = [y1, y2,..,yn]' waypoints y-coord. expressed in NED (m)
%  [x, y]: Vehicle x-y coordinates (m), used to compute y_e
%  h: sampling time (s)
%  undulation_factor: A parameter larger than 0 (typically 1)that controls  
%    the smoothness of the spline's transition between waypoints by scaling
%    the tangents. The undulation_factor directly affects the spline's
%    curvature and reduce undulation, which again minimize the turning rate.
%  Delta_h: Look-ahead distance (m)
%  gamma_h: Positive adaptive gain constant, only used for psi_ref. 
%
% Outputs:
%
%  LOSangle: desired course angle chi_d or desired yaw angle psi_d (rad)
%  y_e:      cross-track error expressed in the path-tangential frame (m)
%  pi_h:     path-tangential angle with respect to NED (rad)
%  closestPoint:   closest point on the current segment to the vehicle
%  closestTangent: tangent at the closest point
%
% Calls: hermiteSpline.m
%
% Reference:
% T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
% Control. 2nd. Edition, Wiley
%
% See also: crosstrackWpt.m, LOSchi.m, ALOSchi.m, ALOSpsi.m
% and SIMotter.m and demoOtterUSVPathFollowingHeadingControl.slx
% for example implementations using an Otter USV.
%
% Author:    Thor I. Fossen
% Date:      2024-04-01
% Revisions: 

persistent beta_hat;  % estimate of the crab angle

if isempty(beta_hat)
    beta_hat = 0;     % initial parameter estimate 
end

if nargin == 5        % course control, no parameter adaptation
    gamma_h = [];
end

n = size(wayPoints, 1);

% Initialize minimum distance to a large value
minDistance = inf;
closestPoint = [0, 0];
closestTangent = [0, 0];
side = 0;

% Add one dummy waypoint, bearing defined by the last two waypoints
R_inf = 10 * max(max(wayPoints));
bearing = atan2(wayPoints(n,2)-wayPoints(n-1,2), ...
    wayPoints(n,1)-wayPoints(end-1,1));
wayPoints = [wayPoints
    wayPoints(n,1) + R_inf * cos(bearing)...
    wayPoints(n,2) + R_inf * sin(bearing)  ];

% Loop through each segment of the Hermite spline
for i = 1:n

    % Define the function to find the closest point on the current
    % segment to the vehicle
    distanceFunc = @(t) norm(hermiteSpline(t, i, wayPoints,...
        undulation_factor) - vehiclePos);

    % Find the t value that minimizes the distance for the current segment
    options = optimset('TolX',1e-6,'Display','off');
    t_min = fminbnd(distanceFunc, 0, 1, options);

    % Calculate the minimum distance for the current segment
    [pointOnSpline, tangentAtPoint] = ...
        hermiteSpline(t_min, i, wayPoints, undulation_factor);
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

if nargin == 5 % LOS guidance law for course control
    
    % Desired course angle: chi_d  
    LOSangle = pi_h - atan( y_e/Delta_h ); 

else % ALOS guidance law for heading control
    
    % Desired heading angle: psi_d  
    LOSangle = pi_h - beta_hat - atan( y_e/Delta_h );

    % Crab angle parameter estimate: beta_hat[k+1]
    beta_hat = beta_hat + ...
        h * gamma_h * Delta_h * y_e / sqrt( Delta_h^2 + y_e^2 );

end

end