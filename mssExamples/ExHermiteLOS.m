% ExHermiteLOS is a example script for computation of the vehicle cross-
% track error, the Line-of-Sight (LOS) course angle and course rate 
% commands when the path isparametrized using a cubic Hermite spline 
% through given waypoints. 
%
% This script visualizes the path, vehicle position, and navigational vectors 
% including the LOS vector, cross-track error, and tangent vectors at
% closest points on the path. 
%
% Calls:      crossTrackErrorHermite.m, hermiteSpline.m     
%
% Author:     Thor I. Fossen
% Date:       2024-02-28
clearvars;

% Waypoint table (xy-coordinates)
wayPoints = [ 0 20
              9 50 
              30 60
              35 80
              55 100];

% Initial vehicle position (xy-coordinates)
vehiclePos = [-10 30];

% LOS parameters
Delta_h = 10;   % look-ahead distance
U = 0.5;        % vehicle speed, should be reduced if omega_chi_d is to large

%% Plot the Hermite spline
figure(1); figure(gcf)

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
table = zeros(N,4);

for step = 1:N  % Loop to update vehicle position and compute cross-track error

    % Update vehicle position 
    vehiclePos = vehiclePos + [step * 0.25, step * 0.2];
    
    % Compute the cross-track error, closest point on the spline, and tangent
    [y_e, pi_h, chi_d, omega_chi_d, closestPoint, closestTangent] = ...
        crosstrackHermiteLOS(wayPoints, vehiclePos, Delta_h, U);

    table(step,:) = [y_e, pi_h, chi_d, omega_chi_d];

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
    quiver(vehiclePos(2), vehiclePos(1), R*sin(chi_d), R*cos(chi_d), 0, ...
        'k', 'LineWidth', 1.5);

end

legend({'Hermite Spline', 'Waypoints', 'Vehicle Position', ...
    'Cross-Track Error', 'Tangent Vector', 'LOS Vector'}, 'Location', 'best');
xlabel('Y (East)'); ylabel('X (North)'); axis equal; grid on; hold off;
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Plot the LOS variables
figure(2); figure(gcf);
subplot(221)
plot(1:N,table(:,1), 'LineWidth', 2), grid
title('Cross-track error y_e^p (m)')
subplot(222)
plot(1:N,rad2deg(table(:,2)), 'LineWidth', 2), grid
title('Path-tangential angle \pi_h (deg)')
subplot(223)
plot(1:N,rad2deg(table(:,3)), 'LineWidth', 2), grid
title('LOS course angle \chi_d (deg)')
subplot(224)
plot(1:N,rad2deg(table(:,4)), 'LineWidth', 2), grid
title(['LOS course rate \omega_{\chi_d} (deg/s) for U = ' num2str(U) ' m/s'])
set(findall(gcf,'type','text'),'FontSize',14)


