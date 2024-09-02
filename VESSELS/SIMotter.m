function SIMotter()
% SIMotter is compatibel with MATLAB and GNU Octave (www.octave.org).
% This script simulates the Otter Uncrewed Surface Vehicle (USV) under 
% various control strategies to handle path following in the presence of 
% ocean currents. This script allows the user to select from several control
% methods and simulatesthe USV's performance using a cubic Hermite spline
% or straight-line paths.
%
% The simulation covers:
%   1. PID heading autopilot without path following.
%   2. Adaptive Line-of-Sight (ALOS) control for path following using
%      straight lines and waypoint switching.
%   3. Integral Line-of-Sight (ILOS) control for path following using
%      straight lines and waypoint switching.
%   4. ALOS control for path following using Hermite spline interpolation.
%
% Dependencies:
%   otter.m                 - Dynamics of the Otter USV
%   refModel.m              - Reference model for autopilot systems
%   ALOSpsi.m               - ALOS guidance algorithm for path following
%   ILOSpsi.m               - ILOS guidance algorithm for path following
%   hermiteSpline.m         - Cubic Hermite spline computations
%   crosstrackHermiteLOS.m  - Cross-track error and LOS guidance law for
%                             cubic Hermite splines
%   LOSobserver.m           - Observer for LOS guidance 
%   controlMethods.m        - Menu for choosing control law. 
%
% Simulink Models:
%   demoOtterUSVPathFollowingHeadingControl.slx
%   demoOtterUSVPathFollowingCourseControl.slx
%
% References:
%   T. I. Fossen and A. P. Aguiar (2024). A Uniform Semiglobal Exponential
%   Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path Following.
%   Automatica, 163, 111556. https://doi.org/10.1016/j.automatica.2024.111556
%
%   T. I. Fossen (2023). An Adaptive Line-of-sight (ALOS) Guidance Law for
%   Path Following of Aircraft and Marine Craft. IEEE Transactions on Control
%   Systems Technology, 31(6), 2887-2894. 
%   https://doi.org/10.1109/TCST.2023.3259819
%
% Author: Thor I. Fossen
% Date: 2021-04-25
% Revisions:
%   2023-10-14: Added heading autopilot and reference model.
%   2024-04-01: Integrated ALOS/ILOS path-following control algorithms for
%               straight-line paths and Hermite splines.
%   2024-04-21: Enhanced compatibility with GNU Octave.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method
%               with RK4.

clearvars;                                  % Clear variables from memory
close all;                                  % Close all figure windows
clear ALOSpsi ILOSpsi crosstrackHermiteLOS  % Clear persistent variables

%% USER INPUTS
h  = 0.05;                       % Sampling time [s]
T_final = 1000;	                 % Final simulation time [s]

% Load condition
mp = 25;                         % Payload mass (kg), maximum value 45 kg
rp = [0.05 0 -0.35]';            % Location of payload (m)

% Ocean current
V_c = 0.3;                       % Ocean current speed (m/s)
beta_c = deg2rad(30);            % Ocean current direction (rad)

% Waypoints
wpt.pos.x = [0 0   150 150 -100 -100 200];
wpt.pos.y = [0 200 200 -50  -50  250 250];

% Add intermediate waypoints along the line segments between for better resolution
wpt = addIntermediateWaypoints(wpt, 2);

% ALOS and ILOS parameters
Delta_h = 10;                    % Look-ahead distance
gamma_h = 0.001;                 % ALOS adaptive gain
kappa = 0.001;                   % ILOS integral gain

% Additional parameter for straight-line path following
R_switch = 5;                    % Radius of switching circle
K_f = 0.3;                       % LOS observer gain

% Initial heading, vehicle points towards next waypoint
psi0 = atan2(wpt.pos.y(2) - wpt.pos.y(1), wpt.pos.x(2) - wpt.pos.x(1));

% Additional parameters for Hermite spline path following
Umax = 2;                        % Maximum speed for Hermite spline LOS
idx_start = 1;                   % Initial index for Hermite spline
[w_path, x_path, y_path, dx, dy, pi_h, pp_x, pp_y, N_horizon] = ...
    hermiteSpline(wpt, Umax, h); % Compute Hermite spline for path following

% Otter USV input matrix
[~,~,M, B_prop] = otter();
Binv = invQR(B_prop);            % Invert input matrix for control allocation

% PID heading autopilot parameters (Nomoto model: M(6,6) = T/K)
T = 1;                           % Nomoto time constant
K = T / M(6,6);                  % Nomoto gain constant

wn = 1.5;                        % Closed-loop natural frequency (rad/s)
zeta = 1.0;                      % Closed-loop relative damping factor (-)

Kp = M(6,6) * wn^2;                     % Proportional gain
Kd = M(6,6) * (2 * zeta * wn - 1/T);    % Derivative gain
Td = Kd / Kp;                           % Derivative time constant
Ti = 10 / wn;                           % Integral time constant

% Reference model parameters
wn_d = 1.0;                      % Natural frequency (rad/s)
zeta_d = 1.0;                    % Relative damping factor (-)
r_max = deg2rad(10.0);           % Maximum turning rate (rad/s)

% Propeller dynamics
T_n = 0.1;                       % Propeller time constant (s)
n = [0 0]';                      % Initial propeller speed, [n_left n_right]'

% Initial states
x = zeros(12,1);                 % x = [u v w p q r xn yn zn phi theta psi]'
x(12) = psi0;                    % Heading angle
z_psi = 0;                       % Integral state for heading control
psi_d = psi0;                    % Desired heading angle
r_d = 0;                         % Desired rate of turn
a_d = 0;                         % Desired acceleration

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% Choose control method and display simulation options
methods = {'PID heading autopilot, no path following',...
           'ALOS path-following control using straight lines and waypoint switching',...
           'ILOS path-following control using straight lines and waypoint switching',...
           'ALOS path-following control using Hermite splines'};
ControlFlag = controlMethod(methods);
displayControlMethod(ControlFlag, R_switch, Delta_h);

%% MAIN LOOP
simdata = zeros(nTimeSteps, 14);    % Preallocate table for simulation data

for i = 1:nTimeSteps

    % Measurements with noise
    r = x(6) + 0.001 * randn;       % Yaw rate 
    xn = x(7) + 0.01 * randn;       % North position
    yn = x(8) + 0.01 * randn;       % East position
    psi = x(12) + 0.001 * randn;    % Yaw angle

    % Guidance and control system
    switch ControlFlag
        case 1  % PID heading autopilot with reference model
            % Reference model, step input adjustments
            psi_ref = psi0;
            if t(i) > 100; psi_ref = deg2rad(0); end
            if t(i) > 500; psi_ref = deg2rad(-90); end

            % Reference model propagation
            [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max, ...
                                         zeta_d, wn_d, h, 1);

        case 2  % ALOS heading autopilot straight-line path following
            psi_ref = ALOSpsi(xn, yn, Delta_h, gamma_h, h, R_switch, wpt);
            [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

        case 3  % ILOS heading autopilot straight-line path following
            psi_ref = ILOSpsi(xn, yn, Delta_h, kappa, h, R_switch, wpt);
            [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

        case 4  % ALOS heading autopilot, cubic Hermite spline interpolation
            [psi_ref, idx_start] = crosstrackHermiteLOS(w_path, x_path, ...
                y_path, dx, dy, pi_h, xn, yn, h, Delta_h, pp_x, pp_y, ...
                idx_start, N_horizon, gamma_h);
            [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);
    end

    % PID heading (yaw moment) autopilot and forward thrust
    tau_X = 100;                              % Constant forward thrust
    tau_N = (T/K) * a_d + (1/K) * r_d - ...   % Yaw moment control
            Kp * (ssa(psi - psi_d) + ...      % Proportional term
            Td * (r - r_d) + (1/Ti) * z_psi); % Derivative and integral terms

    % Control allocation
    u = Binv * [tau_X; tau_N];      % Compute control inputs for propellers
    n_c = sign(u) .* sqrt(abs(u));  % Convert to required propeller speeds

    % Store simulation data
    simdata(i, :) = [x', r_d, psi_d];

    % RK4 method x(k+1)
    x = rk4(@otter, h, x, n,mp,rp,V_c,beta_c);

    % Euler's method
    n = n + h/T_n * (n_c - n);              % Update propeller speeds
    z_psi = z_psi + h * ssa(psi - psi_d);   % Update integral state
end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Get screen dimensions
legendSize = 12;

% Simulation data structure
nu   = simdata(:,1:6); 
eta  = simdata(:,7:12); 
r_d = simdata(:,13);    
psi_d = simdata(:,14);  

% Plot positions
figure(1);
set(gcf,'Position',[0.6*scrSz(3),0.2*scrSz(4), ...
    0.5*scrSz(4),0.6*scrSz(4)],'Visible','off');
hold on;
plot(eta(:,2),eta(:,1),'b');  % Plot vehicle position

% Control method specific plots
if ControlFlag == 1  % Heading autopilot
    legend('Vehicle position');
elseif any(ControlFlag == [2, 3])  % Straight line and circles
    plotStraightLinesAndCircles(wpt, R_switch);
else  % Hermite splines
    plotHermiteSplines(y_path, x_path, wpt);
end

xlabel('East (m)', 'FontSize', 14);
ylabel('North (m)', 'FontSize', 14);
title('North-East Positions (m)', 'FontSize', 14);
axis equal;
grid on;
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','legend'),'FontSize',legendSize);
set(1,'Visible', 'on');  % Show figure

% Plot velocities
figure(2);
if ~isoctave; set(gcf, 'Position', [1, 1, 0.3*scrSz(3), scrSz(4)]); end
subplot(611),plot(t,nu(:,1));
xlabel('Time (s)'),title('Surge velocity (m/s)'),grid on;
subplot(612),plot(t,nu(:,2));
xlabel('Time (s)'),title('Sway velocity (m/s)'),grid on;
subplot(613),plot(t,nu(:,3));
xlabel('Time (s)'),title('Heave velocity (m/s)'),grid on;
subplot(614),plot(t,rad2deg(nu(:,4)));
xlabel('Time (s)'),title('Roll rate (deg/s)'),grid on;
subplot(615),plot(t,rad2deg(nu(:,5)));
xlabel('Time (s)'),title('Pitch rate (deg/s)'),grid on;
subplot(616),plot(t,rad2deg(nu(:,6)),t,rad2deg(r_d));
xlabel('Time (s)'),title('Yaw rate (deg/s)'),grid on;
legend('r','r_d');
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','text'),'FontSize',14);
set(findall(gcf,'type','legend'),'FontSize',legendSize);

% Plot speed, heave position and Euler angles
figure(3);
if ~isoctave
    set(gcf, 'Position', [0.3*scrSz(3), 1, 0.3*scrSz(3), scrSz(4)]); 
end
subplot(511),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2));
title('Speed (m/s)'),grid on;
subplot(512),plot(t,eta(:,3),'linewidth',2);
title('Heave position (m)'),grid on;
subplot(513),plot(t,rad2deg(eta(:,4)));
title('Roll angle (deg)'),grid on;
subplot(514),plot(t,rad2deg(eta(:,5)));
title('Pitch angle (deg)'),grid on;
subplot(515),plot(t,rad2deg(unwrap(eta(:,6))),t,rad2deg(unwrap(psi_d)));
xlabel('Time (s)'),title('Yaw angle (deg)'),grid on;
legend('\psi','\psi_d');
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','text'),'FontSize',14);
set(findall(gcf,'type','legend'),'FontSize',legendSize);

% Display the vehicle data and an image of the vehicle
vehicleData = {...
    'Length', '2.0 m',...
    'Beam', '1.08 m',...
    'Draft (no payload)', '13.4 cm',...
    'Draft (25 kg payload)', '19.5 cm',...    
    'Mass (no payload)', '55.0 kg',...
    'Max speed', '3.0 m/s',...
    'Max pos. propeller speed', '993 RPM',...
    'Max neg. propeller speed', '-972 RPM'};
displayVehicleData('Maritime Robotics Otter USV', vehicleData, 'otter.jpg', 4);

end

%% HELPER FUNCTIONS FOR PLOTTING
function plotStraightLinesAndCircles(wpt, R_switch)
    legendLocation = 'best';
    if isoctave; legendLocation = 'northeast'; end
    
    % Plot straight lines for straight-line path following
    for idx = 1:length(wpt.pos.x)-1
        if idx == 1
            plot([wpt.pos.y(idx), wpt.pos.y(idx+1)], [wpt.pos.x(idx), wpt.pos.x(idx+1)], ...
                 'r--', 'DisplayName', 'Straight-line path');
        else
            plot([wpt.pos.y(idx), wpt.pos.y(idx+1)], [wpt.pos.x(idx), wpt.pos.x(idx+1)], ...
                 'r--', 'HandleVisibility', 'off');
        end
    end

    % Plot circles for waypoint switching
    theta = linspace(0, 2*pi, 100);
    for idx = 1:length(wpt.pos.x)
        xCircle = R_switch * cos(theta) + wpt.pos.x(idx);
        yCircle = R_switch * sin(theta) + wpt.pos.y(idx);
        if idx == 1
            plot(yCircle, xCircle, 'k', 'DisplayName', 'Circle of acceptance');
        else
            plot(yCircle, xCircle, 'k', 'HandleVisibility', 'off');
        end
    end

    % Add legend
    legend('Vehicle position', 'Location', legendLocation);
end

function plotHermiteSplines(y_path, x_path, wpt)
    legendLocation = 'best';
    if isoctave; legendLocation = 'northeast'; end

    % Plot Hermite spline paths for spline path following
    plot(y_path, x_path, 'r');
    plot(wpt.pos.y, wpt.pos.x, 'ko', ...
        'MarkerFaceColor', 'g', 'MarkerSize',10);
    legend('Vehicle position','Hermite spline','Waypoints', ...
        'Location',legendLocation);
end

%% DISPLAY CONTROL METHOD
function displayControlMethod(ControlFlag, R_switch, Delta_h)
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Otter USV');
    switch ControlFlag
        case 1
            disp('PID heading autopilot with reference feedforward');
        case 2
            disp(['ALOS path-following control using straight lines and ' ...
                'waypoint switching']);
            disp(['Circle of acceptance: R = ', num2str(R_switch), ' m']);
            disp(['Look-ahead distance: Delta_h = ', num2str(Delta_h), ' m']);
        case 3
            disp(['ILOS path-following control using straight lines and ' ...
                'waypoint switching']);
            disp(['Circle of acceptance: R =  ', num2str(R_switch), ' m']);
            disp(['Look-ahead distance: Delta_h = ', num2str(Delta_h), ' m']);
        case 4
            disp('ALOS path-following control using Hermite splines');
            disp(['Look-ahead distance: Delta_h = ', num2str(Delta_h), ' m']);
    end
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end
