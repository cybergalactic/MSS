% The script computes a C1 Hermite/MAKIMA spline path through planar
% waypoints given in North-East coordinates. The script also simulates a
% vehicle moving at constant speed and constant course, projects the
% vehicle position onto the sampled path, and visualizes the corresponding
% cross-track error.
%
% The script:
%   1) Generates a Hermite/MAKIMA spline path through the waypoints.
%   2) Computes the path tangent angle pi_h(varpi) and curvature kappa(varpi).
%   3) Simulates a vehicle trajectory in the North-East plane.
%   4) Projects the vehicle position locally onto the path.
%   5) Plots the path, vehicle trajectory, projected points, and the
%      orthogonal cross-track error lines.
%
% Inputs:
%   The script uses internally defined:
%     wpt.pos.x  - Waypoint North coordinates (m)
%     wpt.pos.y  - Waypoint East coordinates (m)
%     Umax       - Maximum/path-following speed used in the spline setup (m/s)
%     h          - Sampling time for vehicle simulation (s)
%
% Outputs:
%   Figure 1:
%     - Hermite spline path
%     - Waypoints
%     - Simulated vehicle trajectory
%     - Projected points on the path
%     - Cross-track error lines
%
%   Figure 2:
%     - Path tangent angle pi_h(varpi)
%     - Path curvature kappa(varpi)
%
% Author:    Thor I. Fossen
% Date:      2026-03-30

clearvars;
clf;

Umax = 1;
h = 0.1;

% Waypoints (North-East coordinates)
wpt.pos.x = [0 150 300 400 200];
wpt.pos.y = [0 120 -50 200 300];

% Compute Hermite spline through the waypoints
[w_path, x_path, y_path, dx_path, dy_path, pi_h, pp_x, pp_y, N_horizon, kappa] = ...
    hermiteSpline(wpt, Umax, h);

% Vehicle trajectory
U = Umax;                 % Vehicle speed (m/s)
dsFwd = 2000;             % Forward search distance
k_old = 1;                % Previous projection index

xN  = 0;                  % Initial North position (m)
yE  = 100;                % Initial East position (m)
chi = deg2rad(20);        % Vehicle course (rad)

Nsim = 5;                 % Number of plotted vehicle points
plotStep = 500;           % Number of simulation steps between plotted points

xVeh   = zeros(Nsim,1);
yVeh   = zeros(Nsim,1);
xProj  = zeros(Nsim,1);
yProj  = zeros(Nsim,1);
x_e    = zeros(Nsim,1);
y_e    = zeros(Nsim,1);
piProj = zeros(Nsim,1);
kProj  = zeros(Nsim,1);

for i = 1:Nsim
    for j = 1:plotStep
        xN = xN + h * U * cos(chi);   % North update
        yE = yE + h * U * sin(chi);   % East update

        % Local projection of vehicle position onto sampled path
        proj = projectToPath(w_path, x_path, y_path, xN, yE, pi_h, k_old, dsFwd);
        k_old = proj.k;               % Update previous projection index
    end

    xVeh(i)   = xN;
    yVeh(i)   = yE;
    xProj(i)  = proj.xp;
    yProj(i)  = proj.yp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; figure(gcf)

subplot(111)
hPath = plot(y_path, x_path, 'LineWidth', 2); hold on;
hWpt  = plot(wpt.pos.y, wpt.pos.x, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% Plot vehicle trajectory (5 points)
hVeh  = plot(yVeh, xVeh, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 6);

% Plot projected points on path
hProj = plot(yProj, xProj, 'kx', 'LineWidth', 1.5, 'MarkerSize', 8);

% Plot cross-track error lines in green
for i = 1:Nsim
    hErr = plot([yVeh(i) yProj(i)], [xVeh(i) xProj(i)], 'g', 'LineWidth', 2);
end

hold off;
axis equal; grid on;
xlabel('East (m)');
ylabel('North (m)');
title('C^1 Hermite spline through waypoints (North-East)');
legend([hPath hWpt hVeh hProj hErr], ...
       {'Path','Waypoints','Vehicle trajectory','Projected points','Cross-track error'}, ...
       'Location','best');

% Path tangent and curvature plots
figure(2); clf; figure(gcf)

subplot(211)
plot(w_path, rad2deg(pi_h), 'LineWidth', 2);
grid on;
xlabel('\varpi');
ylabel('deg');
title('Path tangent \pi_h(\varpi)');

subplot(212)
plot(w_path, kappa, 'LineWidth', 2);
grid on;
xlabel('\varpi');
ylabel('\kappa');
title('Path curvature \kappa(\varpi)');
