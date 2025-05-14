function SIMnpsauv()
% SIMnpsauv is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script simulates the Naval Postgraduate School (NPS) Autonomous
% Underwater Vehicle (AUV), length 5.3 m, under depth and heading control 
% while exposed to ocean currents. It supports control strategies such as
% MIMO PID pole-placement control and Adaptive Line-of-Sight (ALOS) 
% guidance for 3-D path following.
%
% The simulation covers:
%   1. MIMO PID pole placement for pitch and yaw control (autopilots).
%   2. Adaptive Line-of-Sight (ALOS) control for path following using
%      straight lines and waypoint switching.
%
% Dependencies:
%   npsauv.m              - Dynamics of the NPS AUV
%   refModel.m            - Reference model for autopilot systems
%   ALOS.m                - ALOS guidance algorithm for path following
%   LOSobserver.m         - Observer for LOS guidance 
%   ssa.m                 - Angle wrapping
%
% Simulink Models:
%   demoNPSAUV.slx : Simulink model demonstrating PID heading control.
%
% References:
% References:
%   E. M. Coates and T. I. Fossen (2025). Aspherical Amplitude–Phase Formulation 
%       for 3-D Adaptive Line-of-Sight (ALOS) Guidance with USGES Stability
%       Guarantees, Automatica, Submitted. 
%   T. I. Fossen & P. Aguiar (2024). A Uniform Semiglobal Exponential  
%      Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path 
%      Following. Automatica, 163, 111556. 
%      https://doi.org/10.1016/j.automatica.2024.111556
%   A. J. Healey and Lienard, D. (1993). Multivariable Sliding Mode Control 
%     for Autonomous Diving and Steering of Unmanned Underwater Vehicles,
%     IEEE Journal of Ocean Engineering 18(3):327-339.
%
% Author: Thor I. Fossen
% Date: 2024-06-05
% Revisions:
%   2024-07-10: Improved numerical accuracy by replacing Euler's method with RK4.
%   2025-05-13: Crab angle plots based on spherical representation (Coates and 
%               Fossen 2025).

clearvars;                          % Clear all variables from memory
clear ALOS3D;                       % Clear persistent states in controllers
close all;                          % Close all open figure windows

%% USER INPUTS
T_final = 1500;	                    % Final simulation time (s)
h = 0.05;                           % Sampling time (s)

ControlFlag = controlMethod();

% Define waypoints for 3D path following
wpt.pos.x = [0   50  100   0  100, 200  400];
wpt.pos.y = [0  200  600 950 1300 1800 2200];
wpt.pos.z = [0   10  100 200  200  200  150];

% Initialize position and orientation
xn = 0; yn = 0; zn = 0;             % Initial North-East-Down positions (m)
phi = 0; theta = 0;                 % Initial Euler angles (radians)
psi = atan2(wpt.pos.y(2) - wpt.pos.y(1), ...
    wpt.pos.x(2) - wpt.pos.x(1));   % Yaw angle towards next waypoint
U = 1;                              % Initial speed (m/s)

% Initial control and state setup
theta_d = 0; q_d = 0;               % Initial pitch references
psi_d = psi; r_d = 0; a_d = 0;      % Initial yaw references

% Initialize ocean current parameters
Vc = 0.3;                      % Horizontal speed (m/s)
betaVc = deg2rad(20);          % Horizontal direction (radians)
wc = 0.1;                      % Vertical speed (m/s)

% Initialize propeller dynamics
n = 1000;                      % Initial propeller speed (rpm)
n_d = 1300;                    % Desired propeller speed (rpm)
rangeCheck(n_d, 0, 1500);      % Check if within operational limits

% Intitial state vector
x = [U; zeros(5,1); xn; yn; zn; phi; theta; psi; 0; 0; 0; 0; n];

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

%% UNCONSTRAINED CONTROL ALLOCATION
% [tau2 tau3 tau5 tau6] = B_delta * [delta_r, delta_s, delta_bp, delta_bs]
[~, ~, M, B_delta] = npsauv();  % Mass matrix M and input matrix B_delta

% Pseudoinverse (Fossen 2021, Section 11.2.2)
% [delta_r, delta_s, delta_bp, delta_bs] = B_pseudo * [tau5, tau6]
W = diag([5 5 1 1]);   % 5 times more expensive to use delta_r and delta_s
B_pseudo = inv(W) * B_delta' * inv(B_delta * inv(W) * B_delta');

%% CONTROL SYSTEM CONFIGURATION
% Setup for depth and heading control
psi_step = deg2rad(-60);       % Step change in heading angle (rad)
z_step = 30;                   % Step change in depth, max 1000 m
rangeCheck(z_step,0,1000);

% Integral states for autopilots
z_int = 0;                     % Integral state for depth control
theta_int = 0;                 % Integral state for pitch control
psi_int = 0;                   % Integral state for yaw control

% Depth controller (suceessive-loop closure)
z_d = zn;                      % Initial depth target (m)
wn_d_z = 0.02;                 % Natural frequency for depth control
Kp_z = 0.1;                    % Proportional gain for depth
T_z = 100;                     % Time constant for integral action in depth control

% Closed-loop pitch and heading control parameters (rad/s)
zeta_theta = 1.0;              % Damping ratio for pitch control (-)
wn_theta = 1.2;                % Natural frequency for pitch control (rad/s)
zeta_psi = 1.0;                % Damping ratio for yaw control (-)
wn_psi = 0.8;                  % Natural frequency for yaw control (rad/s)

% Heading autopilot reference model parameters
zeta_d_psi = 1.0;              % Damping ratio for yaw control (-)
wn_d_psi = 0.1;                % Natural frequency for yaw control (rad/s)
r_max = deg2rad(10.0);         % Maximum turning rate (rad/s)

% MIMO PID pole-placement algorithm (Algorithm 15.2 in Fossen 2021)
Omega_n = diag([wn_theta wn_psi]);
Zeta = diag([zeta_theta zeta_psi]);
M = diag([M(5,5), M(6,6)]);
Kp = M .* Omega_n.^2;              % Proportional gain
Kd = M .* (2 * Zeta .* Omega_n);   % Derivative gain
Ki = (1/10) * Kp .* Omega_n;       % Integral gain 

%% ALOS PATH-FOLLOWING PARAMETERS
Delta_h = 20;               % horizontal look-ahead distance (m)
Delta_v = 20;               % vertical look-ahead distance (m)
gamma_h = 0.001;            % adaptive gain, horizontal plane
gamma_v = 0.001;            % adaptive gain, vertical plane
M_theta = deg2rad(20);      % maximum value of estimates, alpha_c, beta_c

% Additional parameter for straigh-line path following
R_switch = 5;               % radius of switching circle
K_f = 0.4;                  % LOS observer gain

%% MAIN LOOP
simdata = zeros(nTimeSteps, 29); % Preallocate table for simulation data
ALOSdata = zeros(nTimeSteps, 4); % Preallocate table for ALOS guidance data

for i = 1:nTimeSteps

    % Measurement updates
    u = x(1);                  % Surge velocity (m/s)
    %v = x(2);                 % Sway velocity (m/s)
    %w = x(3);                 % Heave velocity (m/s)
    q = x(5);                  % Pitch rate (rad/s)
    r = x(6);                  % Yaw rate (rad/s)
    xn = x(7);                 % North position (m)
    yn = x(8);                 % East position (m)
    zn = x(9);                 % Down position (m), depth
    %phi = x(10);              % Roll angle (rad) 
    theta = x(11);             % PitchYaw angle (rad)
    psi = x(12);               % Yaw angle (rad)
    %delta_r = x(13);          % Rudder angle (rad)
    %delta_s = x(14);          % Stern plane (rad)
    %delta_bp = x(15);         % Port bow plane (rad)
    %delta_bs  = x(16);        % Starboard bow plane (rad)
    n = x(17);                 % Propeller shaft speed (rpm)

    % Control system updates based on selected mode
    if ControlFlag == 1 % Depth control

        % Depth command, z_ref
        if t(i) > 200
            z_ref = z_step;
        else
            z_ref = 10;
        end

        % LP filtering the depth command
        Uv = sqrt( x(1)^2 + x(3)^2);    % Vertical speed
        if Uv < 1.0                     % Reduce bandwidth at low speed
            wnz = wn_d_z / 2;
        else
            wnz = wn_d_z;
        end
        z_d = exp(-h*wnz) * z_d + (1 - exp(-h*wnz)) * z_ref;

        % Depth autopilot pitch command (succesive-loop closure)
        theta_d = Kp_z * ( (zn - z_d) + (1/T_z) * z_int );     % PI

        % PID heading angle command, psi_ref
        if t(i) > 200
            psi_ref = psi_step;
        else
            psi_ref = deg2rad(0);
        end

        % Third-order reference model for the heading angle
        [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
            zeta_d_psi, wn_d_psi, h, 1);

    else % ALOS path-following

        % ALOS guidance law
        [psi_ref, theta_ref, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
            ALOS3D(xn, yn, zn, Delta_h, Delta_v, gamma_h, gamma_v,...
            M_theta, h, R_switch, wpt);

        % ALOS observer
        [theta_d, q_d] = LOSobserver(theta_d, q_d, theta_ref, h, K_f);
        [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

        ALOSdata(i,:) = [y_e z_e alpha_c_hat beta_c_hat];

    end

    % MIMO PID controller for pitch and roll moments
    tau5 = -Kp(1,1) * ssa( theta - theta_d ) -Kd(1,1) * q ...
        - Ki(1,1) * theta_int;
    tau6 = -Kp(2,2) * ssa( psi - psi_d ) -Kd(2,2) * r ...
        - Ki(2,2) * psi_int;

    % Propeller command (RPM)
    if (n < n_d)
        n = n + 1;
    end

    % Control inputs
    u_com = [ (1/u^2) * B_pseudo * [tau5 tau6]'
        n                                    ];
    max_u = [deg2rad(20); deg2rad(20); deg2rad(20); deg2rad(20); 1500];
    u_com = sat(u_com, max_u);

    % Ocean current random walks with saturation
    betaVc = sat(betaVc + 0.01 * h * randn, deg2rad(100));
    Vc = sat(Vc + 0.05 * h * randn, 1.0);
    wc = sat(wc + 0.01 * h * randn, 0.2);

    % Store simulation data in a table
    simdata(i,:) = [z_d theta_d psi_d r_d u_com' x' Vc betaVc wc];

    % RK methhod (k+1)
    x = rk4(@npsauv, h, x, u_com, Vc, betaVc, wc);  % AUV dynamics 

    % Euler's integration method (k+1)
    z_int = z_int + h * ( zn - z_d );
    theta_int = theta_int + h * ssa( theta - theta_d );
    psi_int = psi_int + h * ssa( psi - psi_d );

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]
legendLocation = 'best'; legendSize = 12;
if isoctave; legendLocation = 'northeast'; end

% simdata(i,:) = [z_d theta_d psi_d r_d u_com' x' Vc betaVc wc];
z_d      = simdata(:,1);
theta_d  = simdata(:,2);
psi_d    = simdata(:,3);
r_d      = simdata(:,4);
u_com    = simdata(:,5:9);
nu       = simdata(:,10:15);
eta      = simdata(:,16:21);
u_actual = simdata(:,22:26);
Vc       = simdata(:,27);
betaVc   = simdata(:,28);
wc       = simdata(:,29);
u = nu(:,1); v = nu(:,2); w = nu(:,3);
p = nu(:,4); q = nu(:,5); r = nu(:,6);


% ALOSdata = [y_e z_e alpha_c_hat beta_c_hat]
y_e = ALOSdata(:,1);
z_e = ALOSdata(:,2);
alpha_c_hat = ALOSdata(:,3);
beta_c_hat = ALOSdata(:,4);

uc = Vc .* cos(betaVc);
vc = Vc .* sin(betaVc);

% Crab angles, AOA and SSA
U = sqrt(u.^2+v.^2+w.^2); % Speed
gamma = asin( (u.*sin(theta)-(v.*sin(phi)+w.*cos(phi)).*cos(theta)) ./ U );
alpha_c = theta - gamma; % Horizontal crab angle
beta_c = atan2(v.*cos(phi)-w.*sin(phi), ...
    u.*cos(theta)+(v.*sin(phi)+w.*cos(phi)).*sin(theta)); % Vertical crab angle
alpha = asin( (w-wc) ./ (u-uc) ); % AOA
beta  = atan2( (v-vc), (u-uc) ); % SSA

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
if ControlFlag == 2; z_d = eta(:,3); end
subplot(411),plot(t,eta(:,3),t,z_d)
xlabel('Time (s)'),title('Down position (m)'),grid
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
if ~isoctave; set(gcf,'Position',[2*scrSz(3)/3,scrSz(4)/2,scrSz(3)/3,scrSz(4)]);end
subplot(511),plot(t,rad2deg(u_actual(:,1)),t,rad2deg(u_com(:,1)))
xlabel('Time (s)'),title('Rudder \delta_r (deg)'),grid
subplot(512),plot(t,rad2deg(u_actual(:,2)),t,rad2deg(u_com(:,2)))
xlabel('Time (s)'),title('Stern-plane \delta_s (deg)'),grid
subplot(513),plot(t,rad2deg(u_actual(:,3)),t,rad2deg(u_com(:,3)))
xlabel('Time (s)'),title('Port bow plane \delta_{bp} (deg)'),grid
subplot(514),plot(t,rad2deg(u_actual(:,4)),t,rad2deg(u_com(:,4)))
xlabel('Time (s)'),title('Starboard bow plane \delta_{bs} (deg)'),grid
subplot(515),plot(t,u_actual(:,5),t,u_com(:,5))
xlabel('Time (s)'),title('Propeller speed n (rpm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

%% Ocean currents and speed
figure(4);
if ~isoctave; set(gcf,'Position',[2*scrSz(3)/3,1,scrSz(3)/4,scrSz(4)/2]);end
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

%% Crab angles, SSA and AOA
if ControlFlag == 2
    figure(5);
    if ~isoctave; set(gcf,'Position',[100,scrSz(4)/2,scrSz(3)/3,scrSz(4)]); end
    subplot(311)
    plot(t,rad2deg(alpha),'g',t,rad2deg(alpha_c),'b',...
        t,rad2deg(alpha_c_hat),'r')
    title('Vertical crab angle and AOA (deg)')
    xlabel('Time (s)')
    grid
    legend('\alpha Angle of attack (AOA)','\alpha_c Vertical crab angle','\alpha_c ALOS estimate','Location',legendLocation)
    subplot(312)
    plot(t,rad2deg(beta),'g',t,rad2deg(beta_c),'b',...
        t,rad2deg(beta_c_hat),'r')
    title('Horizontal crab angle and SSA (deg)')
    xlabel('Time (s)')
    grid
    legend('\beta Sideslip angle (SSA)','\beta_c Horizontal crab angle','\beta_c ALOS estimate','Location',legendLocation)
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
if ControlFlag == 2
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
    ylim([0, 250]);
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
if ControlFlag == 2
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
    'Length', '5.3 m',... 
    'Mass', '5443 kg',...
    'Max speed', '1.9 m/s',...
    'Max propeller speed', '1500 RPM'};
displayVehicleData('Naval Postgraduate School AUV', vehicleData, 'npsauv.png', 9);


end

%% FUNCTIONS
function ControlFlag = controlMethod()

f = figure('Position', [400, 400, 400, 200], ...
    'Name', 'Control Method', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal');

% Add button group for control methods
bg1 = uibuttongroup('Parent', f, ...
    'Position', [0.02 0.6 0.96 0.5], ...
    'Title', 'Control Methods', ...
    'FontSize',14, ...
    'FontWeight','bold');
radio1 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'PID pole-placement control', ...
    'Position', [10 40 500 30], ...
    'Tag', '1');
radio2 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'ALOS guidance law for 3-D path following', ...
    'Position', [10 10 500 30], ...
    'Tag', '2');

% Add OK button to confirm selections
uicontrol('Style', 'pushbutton', ...
    'String', 'OK', ...
    'FontSize', 13, ...
    'Position', [20 50 100 40], ...
    'Callback', @(src, evt) uiresume(f));

uiwait(f); % wait for uiresume to be called on figure handle

% Determine which control method was selected
if get(radio1, 'Value') == 1
    ControlFlag = str2double(get(radio1, 'Tag'));
else
    ControlFlag = str2double(get(radio2, 'Tag'));
end

close(f);  % close the figure after obtaining the selections

disp('-------------------------------------------------------------');
disp('MSS toolbox: NPS AUV');

if (ControlFlag == 1)
    disp('MIMO PID pole placement for pitch and yaw control');
    disp('Successive-loop closure for depth control');
else
    disp('MIMO PID pole placement for pitch and yaw control');
    disp('ALOS guidance law for 3-D path following')
end
disp('-------------------------------------------------------------');
disp('Simulating...');

end
