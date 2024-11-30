function SIMremus100()
% SIMremus100 is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script simulates the Remus 100 Autonomous Underwater Vehicle (AUV) 
% under depth and heading control while exposed to ocean currents. It 
% supports both Euler angle and unit quaternion kinematics representations 
% and features advanced control strategies including PID pole-placement, 
% integral sliding mode control for heading, and Adaptive Line-of-Sight 
% (ALOS) guidance for 3-D path following.
%
% Dependencies:
%   remus100.m            - Dynamics of the Remus 100 AUV
%   integralSMCheading.m  - Integral sliding mode control for heading
%   refModel.m            - Reference model for autopilot systems
%   ALOS.m                - ALOS guidance algorithm for path following
%   LOSobserver.m         - Observer for LOS guidance 
%   lowPassFilter.m       - Low-pass filter
%   q2euler.m, ssa.m      - Utilities for quaternion to Euler conversion 
%                           and angle wrapping
%   displayVehicleData    - Display the vehicle data and an image of the
%                           vehicle
%
% Simulink Models:
%   demoAUVdepthHeadingControl.slx : Depth and heading control.
%
% Reference:
%   T. I. Fossen & P. Aguiar (2024). A Uniform Semiglobal Exponential  
%      Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path 
%      Following. Automatica, 163, 111556. 
%      https://doi.org/10.1016/j.automatica.2024.111556

% Author: Thor I. Fossen
% Date: 2021-06-28
% Revisions:
%   2022-02-01: Autopilots redesign.
%   2022-05-06: Tuning update for remus100 dynamics.
%   2022-05-08: Support for quaternion kinematics.
%   2024-04-02: Simulation options for ALOS path following.
%   2024-04-21: Extended compatibility with GNU Octave.
%   2024-06-07: Included display of vehicle data using displayVehicleData.m.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method
%               with RK4 for the Euler angle representation.

clearvars;                          % Clear all variables from memory
clear integralSMCheading ALOS3D;    % Clear persistent states in controllers
close all;                          % Close all open figure windows

%% USER INPUTS
h = 0.05;                           % Sampling time [s]
T_final = 1400;	                    % Final simulation time [s]

[ControlFlag, KinematicsFlag] = controlMethod();

% Define waypoints for 3D path following
wpt.pos.x = [0  -20 -100   0  200, 200  400];
wpt.pos.y = [0  200  600 950 1300 1800 2200];
wpt.pos.z = [0   10  100 100   50   50   50];

% Initialize position and orientation
xn = 0; yn = 0; zn = 0;             % Initial North-East-Down positions (m)
phi = 0; theta = 0;                 % Initial Euler angles (radians)
psi = atan2(wpt.pos.y(2) - wpt.pos.y(1), ...
    wpt.pos.x(2) - wpt.pos.x(1));   % Yaw angle towards next waypoint
U = 1;                              % Initial speed (m/s)

% Initial control and state setup
theta_d = 0; q_d = 0;               % Initial pitch references
psi_d = psi; r_d = 0; a_d = 0;      % Initial yaw references
xf_z_d = zn;                        % Initial low-pass filter state

% State vector initialization
if KinematicsFlag == 1 
    % Using Euler angles
    x = [U; zeros(5,1); xn; yn; zn; phi; theta; psi];
else
    % Using quaternion representation
    quat = euler2q(phi, theta, psi);
    x = [U; zeros(5,1); xn; yn; zn; quat];
end

% Initialize ocean current parameters
Vc = 0.5;                      % Horizontal speed (m/s)
betaVc = deg2rad(30);          % Horizontal direction (radians)
wc = 0.1;                      % Vertical speed (m/s)

% Initialize propeller dynamics
n = 1000;                      % Initial propeller speed (rpm)
n_d = 1300;                    % Desired propeller speed (rpm)
rangeCheck(n_d, 0, 1525);      % Check if within operational limits

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

%% CONTROL SYSTEM CONFIGURATION
% Setup for depth and heading control
psi_step = deg2rad(-60);       % Step change in heading angle (rad)
z_step = 30;                   % Step change in depth, max 100 m
rangeCheck(z_step,0,100);

% Integral states for autopilots
z_int = 0;                     % Integral state for depth control
theta_int = 0;                 % Integral state for pitch control
psi_int = 0;                   % Integral state for yaw control

% Depth controller (suceessive-loop closure)
z_d = zn;                      % Initial depth target (m)
wn_d_z = 0.02;                 % Natural frequency for depth control
Kp_z = 0.1;                    % Proportional gain for depth
T_z = 100;                     % Time constant for integral action in depth control
Kp_theta = 5.0;                % Proportional gain for pitch control
Kd_theta = 2.0;                % Derivative gain for pitch control
Ki_theta = 0.3;                % Integral gain for pitch control
K_w = 5.0;                     % Feedback gain for heave velocity

% Heading control parameters (using Nomoto model)
K_yaw = 5 / 20;                % Gain, max rate of turn over max. rudder angle
T_yaw = 1;                     % Time constant for yaw dynamics
zeta_d_psi = 1.0;              % Desired damping ratio for yaw control
wn_d_psi = 0.1;                % Natural frequency for yaw control
r_max = deg2rad(5.0);          % Maximum allowable rate of turn (rad/s)

% Heading autopilot (Equation 16.479 in Fossen 2021)
% sigma = r-r_d + 2*lambda*ssa(psi-psi_d) + lambda^2 * integral(ssa(psi-psi_d))
% delta = (T_yaw*r_r_dot + r_r - K_d*sigma - K_sigma*(sigma/phi_b)) / K_yaw
lambda = 0.1;
phi_b = 0.1;                   % Boundary layer thickness

if ControlFlag == 1% PID control parameters
    K_d = 0.5;                 % Derivative gain for PID controller
    K_sigma = 0;               % Gain for SMC (inactive when using PID)
else                        
    % SMC control parameters
    K_d = 0;                   % Derivative gain inactive in SMC mode
    K_sigma = 0.05;            % Sliding mode control gain
end

%% ALOS PATH-FOLLOWING PARAMETERS
Delta_h = 20;               % horizontal look-ahead distance (m)
Delta_v = 20;               % vertical look-ahead distance (m)
gamma_h = 0.001;            % adaptive gain, horizontal plane
gamma_v = 0.001;            % adaptive gain, vertical plane
M_theta = deg2rad(20);      % maximum value of estimates, alpha_c, beta_c

% Additional parameter for straigh-line path following
R_switch = 5;               % radius of switching circle
K_f = 0.5;                  % LOS observer gain

%% MAIN LOOP
simdata = zeros(nTimeSteps, length(x)+10); % Preallocate table for simulation data
ALOSdata = zeros(nTimeSteps, 4); % Preallocate table for ALOS guidance data

for i = 1:nTimeSteps

    % Measurements 
    u = x(1);                  % Surge velocity (m/s)
    v = x(2);                  % Sway velocity (m/s)
    w = x(3);                  % Heave velocity (m/s)
    q = x(5);                  % Pitch rate (rad/s)
    r = x(6);                  % Yaw rate (rad/s)
    xn = x(7);                 % North position (m)
    yn = x(8);                 % East position (m)
    zn = x(9);                 % Down position (m), depth

   if (KinematicsFlag == 1)
         phi = x(10); theta = x(11); psi = x(12); % Euler angles
   else
         [phi,theta,psi] = q2euler(x(10:13)); % Convert quaternion to Euler angles
   end

    % Control system updates based on selected mode
   if ControlFlag == 1 || ControlFlag == 2   % Depth control 

       % Depth command, z_ref
       if t(i) > 200
           z_ref = z_step;
       else
           z_ref = 10;
       end

       % LP filtering the depth command
       Uv = sqrt(x(1)^2+x(3)^2);    % Vertical speed
       if Uv < 1.0                  % Reduce bandwidth at low speed
           wnz = wn_d_z / 2;
       else
           wnz = wn_d_z;
       end
       [xf_z_d, z_d] = lowPassFilter(xf_z_d, z_ref, wnz, h);

       % Depth autopilot using the stern planes (succesive-loop closure)
       theta_d = Kp_z * ( (zn - z_d) + (1/T_z) * z_int );     % PI
       delta_s = -Kp_theta * ssa( theta - theta_d )...        % PID
           - Kd_theta * q - Ki_theta * theta_int - K_w * w;

       % PID heading angle command, psi_ref
       if t(i) > 200
           psi_ref = psi_step;
       else
           psi_ref = deg2rad(0);
       end

       % Heading autopilot using the tail rudder
       delta_r = integralSMCheading(psi, r, psi_d, r_d, a_d, ...
           K_d, K_sigma, lambda, phi_b, K_yaw, T_yaw, h);

       % Third-order reference model for the heading angle
       [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
           zeta_d_psi, wn_d_psi, h, 1);

   else % ControlFlag == 3, ALOS path-following

       % Heading autopilot using the tail rudder (integral SMC)
       delta_r = integralSMCheading(psi, r, psi_d, r_d, a_d, ...
           K_d, K_sigma, 1, phi_b, K_yaw, T_yaw, h);

       % Depth autopilot using the stern planes (PID)
       delta_s = -Kp_theta * ssa( theta - theta_d )...
           - Kd_theta * q - Ki_theta * theta_int - K_w * w;

       % ALOS guidance law
       [psi_ref, theta_ref, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
           ALOS3D(xn, yn, zn, Delta_h, Delta_v, gamma_h, gamma_v,...
           M_theta, h, R_switch, wpt);

       % ALOS observer
       [theta_d, q_d] = LOSobserver(theta_d, q_d, theta_ref, h, K_f);
       [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);
       if abs(r_d) > r_max, r_d = sign(r_d) * r_max; end

       % Ocean current dynamics
       if t(i) > 800
           Vc_d = 0.65;
           w_V = 0.1;
           Vc = exp(-h*w_V) * Vc + (1 - exp(-h*w_V)) * Vc_d;
       else
           Vc = 0.5;
       end

       if t(i) > 500
           betaVc_d = deg2rad(160);
           w_beta = 0.1;
           betaVc = exp(-h*w_beta) * betaVc + (1 - exp(-h*w_beta)) * betaVc_d;
       else
           betaVc = deg2rad(150);
       end

       betaVc = betaVc + (pi/180) * randn / 20;
       Vc = Vc + 0.002 * randn;

       ALOSdata(i,:) = [y_e z_e alpha_c_hat beta_c_hat];

   end

   % Propeller control (rpm)
   if (n < n_d)
       n = n + 1;
   end

   % Amplitude saturation of the control signals
   n_max = 1525;                                % maximum propeller RPM
   max_ui = [deg2rad(15) deg2rad(15) n_max]';   % deg, deg, RPM

   if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
   if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
   if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end

   ui = [delta_r -delta_s n]';                % Commanded control inputs

   % Store simulation data in a table
   simdata(i,:) = [z_d theta_d psi_d r_d Vc betaVc wc ui' x'];

   if (KinematicsFlag == 1)
       % Euler angles x = [ u v w p q r x y z phi theta psi ]'
       x = rk4(@remus100, h, x, ui, Vc, betaVc, wc);  % RK4 method x(k+1)
   else
       % Unit quaternions x = [ u v w p q r x y z eta eps1 eps2 eps3 ]'  
       xdot = remus100(x, ui, Vc, betaVc, wc);
       quat = x(10:13);                            % Unit quaternion
       x(1:6) = x(1:6) + h * xdot(1:6);            % Forward Euler
       x(7:9) = x(7:9) + h * Rquat(quat) * x(1:3); % Backward Euler
       quat = expm(Tquat(x(4:6)) * h) * quat;  % Exact quat. discretization
       x(10:13) = quat / norm(quat);           % Normalization
   end

   % Euler's integration method (k+1)
   z_int = z_int + h * ( zn - z_d );
   theta_int = theta_int + h * ssa( theta - theta_d );
   psi_int = psi_int + h * ssa( psi - psi_d );

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]
legendLocation = 'best'; legendSize = 12;
if isoctave; legendLocation = 'northeast'; end

% simdata = [z_d theta_d psi_d r_d Vc betaVc wc ui' x']
z_d     = simdata(:,1);
theta_d = simdata(:,2);
psi_d   = simdata(:,3);
r_d     = simdata(:,4);
Vc      = simdata(:,5);
betaVc  = simdata(:,6);
wc      = simdata(:,7);
u       = simdata(:,8:10);
nu      = simdata(:,11:16);

if (KinematicsFlag == 1) % Euler angle representation
    eta = simdata(:,17:22);
else % Transform the unit quaternions to Euler angles
    quaternion = simdata(:,20:23);
    for i = 1:length(t)
        [phi(i,1),theta(i,1),psi(i,1)] = q2euler(quaternion(i,:));
    end
    eta = [simdata(:,17:19) phi theta psi];
end

% ALOSdata = [y_e z_e alpha_c_hat beta_c_hat]
y_e = ALOSdata(:,1);
z_e = ALOSdata(:,2);
alpha_c_hat = ALOSdata(:,3);
beta_c_hat = ALOSdata(:,4);

uc = Vc .* cos(betaVc);
vc = Vc .* sin(betaVc);
alpha_c = atan( (nu(:,2).*sin(eta(:,4))+nu(:,3).*cos(eta(:,4))) ./ nu(:,1) );
Uv = nu(:,1) .* sqrt( 1 + tan(alpha_c).^2 );
beta_c = atan( ( nu(:,2).*cos(eta(:,4))-nu(:,3).*sin(eta(:,4)) ) ./ ...
    ( Uv .* cos(eta(:,5)-alpha_c) ) );
alpha = atan2( (nu(:,3)-wc), (nu(:,1)-uc) );
beta  = atan2( (nu(:,2)-vc), (nu(:,1)-uc) );

%% Generalized velocity
figure(1);
if ~isoctave; set(gcf,'Position',[1, 1, scrSz(3)/3, scrSz(4)]); end
subplot(611),plot(t,nu(:,1))
xlabel('Time (s)'),title('Surge velocity (m/s)'),grid
subplot(612),plot(t,nu(:,2))
xlabel('Time (s)'),title('Sway velocity (m/s)'),grid
subplot(613),plot(t,nu(:,3))
xlabel('Time (s)'),title('Heave velocity (m/s)'),grid
subplot(614),plot(t,(180/pi)*nu(:,4))
xlabel('Time (s)'),title('Roll rate (deg/s)'),grid
subplot(615),plot(t,(180/pi)*nu(:,5))
xlabel('Time (s)'),title('Pitch rate (deg/s)'),grid
subplot(616),plot(t,(180/pi)*nu(:,6))
xlabel('Time (s)'),title('Yaw rate (deg/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Heave position and Euler angles
figure(2);
if ~isoctave; set(gcf,'Position',[scrSz(3)/3, 1, scrSz(3)/3, scrSz(4)]); end
if ControlFlag == 3; z_d = eta(:,3); end
subplot(411),plot(t,eta(:,3),t,z_d)
xlabel('Time (s)'),title('Heave position (m)'),grid
legend('True','Desired')
subplot(412),plot(t,rad2deg(eta(:,4)))
xlabel('Time (s)'),title('Roll angle (deg)'),grid
subplot(413),plot(t,rad2deg(eta(:,5)),t,rad2deg(theta_d))
xlabel('Time (s)'),title('Pitch angle (deg)'),grid
legend('True','Desired')
subplot(414),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d))
xlabel('Time (s)'),title('Yaw angle (deg)'),grid
legend('True','Desired')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Control signals
figure(3);
if ~isoctave; set(gcf,'Position',[2*scrSz(3)/3,scrSz(4)/2,scrSz(3)/3,scrSz(4)/2]);end
subplot(311),plot(t,rad2deg(u(:,1)))
xlabel('Time (s)'),title('Rudder command \delta_r (deg)'),grid
subplot(312),plot(t,rad2deg(u(:,2)))
xlabel('Time (s)'),title('Stern-plane command \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3))
xlabel('Time (s)'),title('Propeller speed command n (rpm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

%% Ocean currents and speed
figure(4);
if ~isoctave; set(gcf,'Position',[2*scrSz(3)/3,1,scrSz(3)/3,scrSz(4)/2]);end
subplot(311),plot(t,sqrt(nu(:,1).^2+nu(:,2).^2),t,Vc)
xlabel('Time (s)'),grid
legend('Vehicle horizontal speed (m/s)','Ocean current horizontal speed (m/s)',...
    'Location',legendLocation)
subplot(312),plot(t,nu(:,3),t,wc)
xlabel('Time (s)'),grid
legend('Vehicle heave velocity (m/s)','Ocean current heave velcoity (m/s)',...
    'Location',legendLocation)
subplot(313),plot(t,rad2deg(betaVc),'r')
xlabel('Time (s)'),grid
legend('Ocean current direction (deg)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Sideslip and angle of attack
if ControlFlag == 3
    figure(5);
    if ~isoctave; set(gcf,'Position',[100,scrSz(4)/2,scrSz(3)/3,scrSz(4)]); end
    subplot(311)
    plot(t,rad2deg(alpha),'g',t,rad2deg(alpha_c),'b',...
        t,rad2deg(alpha_c_hat),'r')
    title('Angle of attack (deg)')
    xlabel('Time (s)')
    grid
    legend('\alpha','\alpha_c','\alpha_c estimate','Location',legendLocation)
    subplot(312)
    plot(t,rad2deg(beta),'g',t,rad2deg(beta_c),'b',...
        t,rad2deg(beta_c_hat),'r')
    title('Sideslip angle (deg)')
    xlabel('Time (s)')
    grid
    legend('\beta','\beta_c','\beta_c estimate','Location',legendLocation)
    subplot(313)
    plot(t,y_e,t,z_e)
    title('Tracking errors (m)'),grid
    xlabel('Time (s)')
    legend('Cross-track error y_e^p','Vertical-track error z_e^p')
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',legendSize)
end

%% 2-D position plots with waypoints
if ControlFlag == 3
    figure(6);
    if ~isoctave;set(gcf,'Position',[300,200,scrSz(3)/3,scrSz(4)/2]);end
    subplot(211);
    plot(eta(:,2), eta(:,1));
    hold on;
    plot(wpt.pos.y, wpt.pos.x, 'rx', 'MarkerSize', 10);
    hold off;
    xlabel('East');
    ylabel('North');
    title('North-East plot (m)');
    xlim([0, 2500]);
    axis('equal');
    grid on;
    subplot(212);
    plot(eta(:,2), eta(:,3));
    hold on;
    plot(wpt.pos.y, wpt.pos.z, 'rx', 'MarkerSize', 10);
    hold off;
    xlim([0, 2500]);
    ylim([0, 150]);
    xlabel('East');
    ylabel('Down');
    title('Down-East plot (m)');
    grid on;
    legend('Actual path', 'Waypoints', 'Location', legendLocation);
    set(findall(gcf, 'type', 'line'), 'LineWidth', 2);
    set(findall(gcf, 'type', 'text'), 'FontSize', 14);
    set(findall(gcf, 'type', 'legend'), 'FontSize', legendSize);
end

%% 3-D position plot with waypoints
if ControlFlag == 3
    figure(7);
    plot3(eta(:,2),eta(:,1),eta(:,3))
    hold on;
    plot3(wpt.pos.y, wpt.pos.x, wpt.pos.z, 'ro', 'MarkerSize', 15);
    hold off
    title('North-East-Down plot (m)')
    xlabel('East'); ylabel('North'); zlabel('Down');
    legend('Actual path','Waypoints','Location',legendLocation),grid
    set(gca, 'ZDir', 'reverse');
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',legendSize)
    view(-25, 30);  % view(AZ,EL)
end

% Display the vehicle data and an image of the vehicle
vehicleData = {...
    'Length', '1.6 m',...
    'Diameter', '19 cm',...
    'Mass', '31.9 kg',...
    'Max speed', '2.5 m/s',...
    'Max propeller speed', '1525 RPM'};
displayVehicleData('Remus100 AUV', vehicleData, 'remus100.jpg', 8);

end

%% FUNCTIONS
function [ControlFlag, KinematicsFlag] = controlMethod()

f = figure('Position', [400, 400, 400, 400], ...
    'Name', 'Control Method and Kinematic Representation', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal');

% Add button group for control methods
bg1 = uibuttongroup('Parent', f, ...
    'Position', [0.02 0.6 0.96 0.3], ...
    'Title', 'Control Methods', ...
    'FontSize',14, ...
    'FontWeight','bold');
radio1 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, ...
    'String', 'PID pole-placement control', ...
    'Position', [10 70 500 30], ...
    'Tag', '1');
radio2 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'Integral sliding mode control (SMC)', ...
    'Position', [10 40 500 30], ...
    'Tag', '2');
radio3 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'ALOS guidance law for 3-D path following', ...
    'Position', [10 10 500 30], ...
    'Tag', '3');

% Add button group for kinematics options
bg2 = uibuttongroup('Parent', f, ...
    'Position', [0.02 0.3 0.96 0.25], ...
    'Title', 'Kinematics Representation', ...
    'FontSize',14, ...
    'FontWeight','bold');
radio4 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 13, 'String', ...
    'Euler angles', ...
    'Position', [10 45 500 30], ...
    'Tag', '1', ...
    'Value', 0);
radio5 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 13, 'String', ...
    'Unit quaternions', ...
    'Position', [10 15 500 30], ...
    'Tag', '2', ...
    'Value', 1);

% Add OK button to confirm selections
uicontrol('Style', 'pushbutton', ...
    'String', 'OK', ...
    'FontSize', 13, ...
    'Position', [20 50 100 40], ...
    'Callback', @(src, evt) uiresume(f));

uiwait(f); % Wait for uiresume to be called on figure handle

% Determine which control method was selected
if get(radio1, 'Value') == 1
    ControlFlag = str2double(get(radio1, 'Tag'));
elseif get(radio2, 'Value') == 1
    ControlFlag = str2double(get(radio2, 'Tag'));
else
    ControlFlag = str2double(get(radio3, 'Tag'));
end

% Determine which kinematics option was selected
if get(radio4, 'Value') == 1
    KinematicsFlag = str2double(get(radio4, 'Tag'));
else
    KinematicsFlag = str2double(get(radio5, 'Tag'));
end

close(f);  % Close the figure after obtaining the selections

disp('-------------------------------------------------------------');
disp('MSS toolbox: Remus 100 AUV');

if (KinematicsFlag == 1)
    disp('Euler angle representation (12 states)');
else
    disp('Unit quaternion representation (13 states)');
end

if (ControlFlag == 1)
    disp('Heading autopilot: PID pole-placement control');
    disp('Depth autopilot:   Successive-loop closure');
elseif (ControlFlag == 2)
    disp('Heading autopilot: Integral sliding mode control (SMC)');
    disp('Depth autopilot:   Successive-loop closure');
else
    disp('Heading autopilot: Integral sliding mode control (SMC)');
    disp('Depth autopilot:   Successive-loop closure');
    disp('Path-following:    ALOS guidance law for 3-D path following')
end
disp('-------------------------------------------------------------');
disp('Simulating...');

end
