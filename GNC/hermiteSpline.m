function [P, T] = hermiteSpline(t, segmentIndex, wayPoints)
% [P, T] = hermiteSpline(t, segmentIndex, wayPoints) computes the cubic 
% Hermite spline P and the tangents T to the spline for n > 3 waypoints.
% The tangents at the waypoints are computed using finite differences, 
% ensuring a smooth and continuous spline that is customizable through the 
% choice of waypoints.
%
% Inputs:    
%  t: parameter (or a vector of parameters) that represents the normalized 
%    position(s) within a segment of the spline. It varies from 0 to 1, where 
%    0 corresponds to the start of the segment and 1 corresponds to the end 
%    the segment. t is used to evaluate the spline and its derivative at 
%    specific points along the segment.
%  segmentIndex: integer that specifies the index of the segment for which 
%    the spline and its tangents are to be computed. The segments are 
%    defined by the waypoints, with each segment spanning from one waypoint 
%    to the next. For n waypoints, there are n-1 segments, and segmentIndex 
%    should be in the range of 1 to n-1.
%  wayPoints: matrix where each row represents a waypoint in 2-D space, 
%    with the first column corresponding to the x-coordinate and the second 
%    column to the y-coordinate. These waypoints define the path along which
%    the Hermite spline is interpolated. The spline will be smooth and 
%    continuous, passing through each waypoint and oriented according to 
%    the computed tangents.
%
% Outputs:  
%  P: matrix or vector (depending on the size of t) that represents the 
%    points on the cubic Hermite spline at the parameter(s) t for the 
%    specified segment. Each row of P corresponds to a point in 2-D space, 
%    with the first column representing the x-coordinate and the second 
%    column the y-coordinate. These points are computed based on the 
%    Hermite basis functions and the positions and tangents of the waypoints 
%    defining the segment.  
%  T: matrix or vector that represents the tangents to the cubic Hermite 
%    spline at the parameter(s) t for the specified segment. Like P, each 
%    row of T corresponds to a tangent vector in 2-D space, with the first
%    column representing the x-component and the second column the
%    y-component. These tangents are derived from the derivatives of the
%    Hermite basis functions and provide information about the direction 
%    and speed of the spline at each point P.
%
% See: [y_e, pi_h, chi, U_max, closestPoint, closestTangent] = ...
%    crosstrackHermiteLOS(wayPoints, vehiclePos, Delta_h, omega_chi_max)
%
% Author:    Thor I. Fossen
% Date:      2024-02-28
% Revisions: 2012-03-29 - Added a tangent scale factor to reduce undulation

n = size(wayPoints, 1);     % Number of waypoints
tangents = zeros(n, 2);     % Initialize tangents array

% Find the maximum absolute value among the waypoints to determine scaling
maxVal = max(abs(wayPoints), [], 'all');

% By scaling the tangents, the influence of each waypoint's direction on 
% the spline's curvature is reduced, which can help in achieving a 
% smoother transition between waypoints, thus reducing undulation.
scale_factor = 1 / (1 + log10(maxVal + 1));

% Compute and scale tangents using finite differences
for i = 1:n
    if i == 1
        tangents(i, :) = (wayPoints(i+1, :) - wayPoints(i, :)) * scale_factor;
    elseif i == n
        tangents(i, :) = (wayPoints(i, :) - wayPoints(i-1, :)) * scale_factor;
    else
        tangents(i, :) = ((wayPoints(i+1, :) - wayPoints(i-1, :)) / 2) * scale_factor;
    end
end

% Hermite spline interpolation for the given segment
P0 = wayPoints(segmentIndex, :);
P1 = wayPoints(segmentIndex + 1, :);
T0 = tangents(segmentIndex, :);
T1 = tangents(segmentIndex + 1, :);

% Hermite basis functions
h00 = 2*t.^3 - 3*t.^2 + 1;
h10 = t.^3 - 2*t.^2 + t;
h01 = -2*t.^3 + 3*t.^2;
h11 = t.^3 - t.^2;

% Compute the spline points
P = h00.*P0 + h10.*T0 + h01.*P1 + h11.*T1;

% Compute the derivatives for the tangent vectors, if needed
dh00 = 6*t.^2 - 6*t;
dh10 = 3*t.^2 - 4*t + 1;
dh01 = -6*t.^2 + 6*t;
dh11 = 3*t.^2 - 2*t;

T = dh00.*P0 + dh10.*T0 + dh01.*P1 + dh11.*T1;

end
