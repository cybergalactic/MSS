% ExHermiteLOS is a sample script for calculating the cross-track error 
% of a vehicle during 2-D path following. It computes the Line-of-Sight (LOS) 
% course angle and course rate commands for a path that is represented by 
% a cubic Hermite spline through specified waypoints.
%
% This script visualizes the 2-D path, vehicle position, and navigational 
% vectors including the LOS vector, cross-track error, and tangent vectors
% at closest points on the path. 
%
% Calls:      crossTrackErrorHermite.m
%             hermiteSpline.m     
%
% Example:    See 'SIMotter.m' for an example with adaptive line-of-sight
%             (ALOS) path-following control using Hermite splines
%
% Author:     Thor I. Fossen
% Date:       2024-04-01

clearvars;

% Table of waypoints (xy-coordinates)
wpt.pos.x = [0, 9, 30, 35, 55]';
wpt.pos.y = [20, 50, 60, 80, 100]';
wayPoints = [wpt.pos.x wpt.pos.y];

% Initial vehicle position (xy-coordinates)
x = -10;
y = 30;
vehiclePos = [x y];

% LOS parameters
Delta_h = 10;               % look-ahead distance
undulation_factor = 0.5;    % smoothness of the spline between waypoints

%% Plot the Hermite spline
figure(1); clf; figure(gcf)

t = linspace(0, 1, 200);
for i = 1:size(wayPoints, 1)-1
    P = zeros(length(t), 2);
    for j = 1:length(t)
        [P(j, :), ~] = hermiteSpline(t(j), i, wayPoints); 
    end
    if i == 1
       plot(P(:, 2), P(:, 1), 'b-', 'LineWidth', 3); 
       hold on;
    else
       plot(P(:, 2), P(:, 1), 'b-', 'LineWidth', 3,'HandleVisibility', 'off'); 
    end
end

% Plot the waypoints
plot(wayPoints(:, 2), wayPoints(:, 1), 'ko', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 15);

%% MAIN LOOP
N = 25;
table = zeros(N,3);

for h = 1:N  % Loop to update vehicle position and cross-track error

    % Update dummy vehicle position 
    vehiclePos = vehiclePos + [h * 0.25, h * 0.2];
    
    % Compute the cross-track error, closest point on the spline, and tangent
    [chi_ref, y_e, pi_h, closestPoint, closestTangent] = ...
      crosstrackHermiteLOS(wayPoints,vehiclePos,h,undulation_factor,Delta_h);
    % Store data 
    table(h,:) = [y_e, pi_h, chi_ref];

    % Plot the vehicle North-East position
    plot(vehiclePos(2), vehiclePos(1),'rs','MarkerFaceColor','r', ...
        'MarkerSize', 10, 'LineWidth',3);

    % Plot the cross-track error line
    plot([vehiclePos(2), closestPoint(2)], ...
        [vehiclePos(1), closestPoint(1)], 'r-', 'LineWidth', 1);

    % Plot the tangent vector at the closest point
    normalizedTangent = closestTangent / norm(closestTangent);
    scaledTangent = normalizedTangent * Delta_h;
    quiver(closestPoint(2), closestPoint(1), ...
        scaledTangent(2), scaledTangent(1), 0, 'c-', 'LineWidth', 1);

    % Plot the LOS vector
    R = sqrt(Delta_h^2 + y_e^2);
    quiver(vehiclePos(2), vehiclePos(1), R*sin(chi_ref), R*cos(chi_ref), 0, ...
        'k', 'LineWidth', 1.5);

end

legend({'Hermite Spline', 'Waypoints', 'Vehicle Position', ...
    'Cross-Track Error', 'Tangent Vector', 'LOS Vector'}, 'Location', 'best');
xlabel('Y (East)'); ylabel('X (North)'); axis equal; grid on; hold off;
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Plot the LOS variables stored in table(:,1:4)
figure(2); clf; figure(gcf);
subplot(221)
plot(1:N,table(:,1), 'LineWidth', 2), grid
title('Cross-track error y_e^p (m)')
subplot(222)
plot(1:N,rad2deg(table(:,2)), 'LineWidth', 2), grid
title('Path-tangential angle \pi_h (deg)')
subplot(212)
plot(1:N,rad2deg(table(:,3)), 'LineWidth', 2), grid
title('LOS course angle \chi_{ref} (deg)')



