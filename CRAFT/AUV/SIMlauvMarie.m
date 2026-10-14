function SIMlauvMarie()
% SIMlauvMarie simulates LAUV Marie at NTNU with depth/heading or ALOS control.
% The 2.15 m LAUV is simulated under depth and heading control while exposed
% to ocean currents and directional waves. It supports MIMO PID pole-placement
% control and Adaptive Line-of-Sight (ALOS) guidance for 3-D path following.
% Wave-induced forces are attenuated with vehicle depth and are typically
% negligible below about 20 m for the default sea state. The attenuation is
% frequency dependent, so very long waves can penetrate more deeply.
% The vehicle dynamics use a unified seakeeping-maneuvering model based on
% MSS-Capytaine hydrodynamic data. The power-based constant hydrodynamic model
% is augmented with nonlinear damping, cross-flow drag, and lift-drag terms.
%
% The simulation covers:
%   1. MIMO PID pole placement for pitch and yaw control (autopilots).
%   2. Adaptive Line-of-Sight (ALOS) control for path following using
%      straight lines and waypoint switching.
%
% Dependencies:
%   lauvMarie.m           - Dynamics of the LAUV Marie
%   refModel.m            - Reference model for autopilot systems
%   ALOS3D.m              - ALOS guidance algorithm for path following
%   LOSobserver.m         - Observer for LOS guidance
%   ssa.m                 - Angle wrapping
%   displayVehicleData.m  - Display vehicle particulars and photograph
%
% References:
%   E. M. Coates and T. I. Fossen (2026). A Spherical Amplitude–Phase Formulation
%       for 3-D Adaptive Line-of-Sight (ALOS) Guidance with USGES Stability
%       Guarantees, Automatica 194, 113149S.
%   T. I. Fossen & P. Aguiar (2024). A Uniform Semiglobal Exponential
%      Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path
%      Following. Automatica, 163, 111556.
%      https://doi.org/10.1016/j.automatica.2024.111556
%
% Author: Thor I. Fossen
% Date: 2026-10-14
% Revisions:

clear ALOS3D lauvMarie;             % Clear persistent model/controller states

% ==============================================================================
% SIMULATOR CONFIGURATION
% ==============================================================================
h = 0.05;                           % Sampling time (s)
T_final = 1000;	                    % Final simulation time (s)

[ControlFlag, Hs, w0_spectrum, beta_wave, Vc, betaVc, wc] = controlMethod();

% ==============================================================================
%% USER INPUTS
% ==============================================================================
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

% Initial 12-state vector [nu; eta]
x = [U; zeros(5,1); xn; yn; zn; phi; theta; psi];

% ==============================================================================
%% WAVE SPECTRUM AND FORCE RAO INITIALIZATIONS
% ==============================================================================
environment.Hs = Hs;
environment.w0 = w0_spectrum;
environment.beta_wave = beta_wave; % Mean wave direction in NED (rad)
environment.spectrumType = 'JONSWAP';
environment.spreadingFlag = true;
environment.numFreqIntervals = 100;
environment.numDirections = 24;
RAO_update_period = 0.1;
rng(1); % Reproducible stochastic wave realization

matFile = resolveHydroVesselFile('LAUV_marie.mat');
auvData = load(matFile, 'vessel');
vessel = auvData.vessel;

[environment, waveModel] = waveInitialization(vessel, environment);

% ==============================================================================
%% CONTROL SYSTEM CONFIGURATION
% ==============================================================================
time_constants_1236 = [50 5 5 5];
delta_zeta_45 = [0.2 0.2];

vessel = computeManeuveringModel(vessel, w0_spectrum, ...
    time_constants_1236, delta_zeta_45, 0);

% Setup for depth and heading control
psi_step = deg2rad(-40);       % Step change in heading angle (rad)
z_step = 30;                   % Step change in depth, max 1000 m
rangeCheck(z_step,0,1000);

% Integral states for autopilots
z_int = 0;                     % Integral state for depth control
theta_int = 0;                 % Integral state for pitch control
psi_int = 0;                   % Integral state for yaw control

% Depth controller (successive-loop closure)
z_d = zn;                      % Initial depth target (m)
wn_d_z = 0.02;                 % Natural frequency for depth control
Kp_z = 0.1;                    % Proportional gain for depth
T_z = 100;                     % Time constant for integral action in depth control
k_grad = 0.1;                  % Gain for computation of theta_d

% Closed-loop pitch and heading control parameters (rad/s)
zeta_theta = 0.8;              % Damping ratio for pitch control (-)
wn_theta = 0.6;                % Natural frequency for pitch control (rad/s)
zeta_psi = 0.8;                % Damping ratio for yaw control (-)
wn_psi = 1.0;                  % Natural frequency for yaw control (rad/s)

% Heading autopilot reference model parameters
zeta_d_psi = 1.0;              % Damping ratio for yaw control (-)
wn_d_psi = 0.05;               % Natural frequency for yaw control (rad/s)
r_max = deg2rad(10.0);         % Maximum turning rate (rad/s)

% MIMO PID pole-placement algorithm (Algorithm 15.2 in Fossen 2027)
Omega_n = diag([wn_theta wn_psi]);
Zeta = diag([zeta_theta zeta_psi]);
M_control = diag([vessel.M(5,5), vessel.M(6,6)]);
Kp = M_control .* Omega_n.^2;              % Proportional gain
Kd = M_control .* (2 * Zeta .* Omega_n);   % Derivative gain
Ki = (1/10) * Kp .* Omega_n;       % Integral gain

%% ALOS PATH-FOLLOWING PARAMETERS
Delta_h = 20;               % horizontal look-ahead distance (m)
Delta_v = 20;               % vertical look-ahead distance (m)
gamma_h = 0.001;            % adaptive gain, horizontal plane
gamma_v = 0.001;            % adaptive gain, vertical plane
M_theta = deg2rad(20);      % maximum value of estimates, alpha_c, beta_c

% Additional parameters for straight-line path following
R_switch = 5;               % radius of switching circle
K_f = 0.4;                  % LOS observer gain

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final
nextRAOtime = 0;                % Next RAO update time
nTimeSteps = length(t);         % Number of time steps

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
simdata = zeros(nTimeSteps, 32); % [references, tau, x, current, waves]
ALOSdata = zeros(nTimeSteps, 4); % Preallocate table for ALOS guidance data
tau_wave1 = zeros(6,1);
waveElevation = 0;

for i = 1:nTimeSteps

    % Measurement updates
    u = x(1);                  % Surge velocity (m/s)
    v = x(2);                  % Sway velocity (m/s)
    w = x(3);                  % Heave velocity (m/s)
    q = x(5);                  % Pitch rate (rad/s)
    r = x(6);                  % Yaw rate (rad/s)
    xn = x(7);                 % North position (m)
    yn = x(8);                 % East position (m)
    zn = x(9);                 % Down position (m), depth
    %phi = x(10);              % Roll angle (rad)
    theta = x(11);             % Pitch angle (rad)
    psi = x(12);               % Yaw angle (rad)
    % Control system updates based on selected mode
    if ControlFlag == 1 % Depth control

        % Depth command, z_ref
        if t(i) > 200
            z_ref = z_step;
        else
            z_ref = 1;
        end

        % Depth kinematics:
        %   z_dot = -u * sin(theta) + w * cos(theta) := alpha_z + sigma
        % Virtual control:
        %   alpha_z = z_d_dot + Kp_z * (e_z + (1/T_z) * z_int)
        %
        % Outer-loop depth control residual sigma = 0 implies e_z = 0
        sigma = -u * sin(theta) + w * cos(theta) ...
            + Kp_z * ((zn - z_d) + (1/T_z) * z_int);

        % Gradient-descent update using J = 1/2 * sigma^2 to drive sigma to 0:
        %   dJ/d(theta) = - u * cos(theta) - w * sin(theta)
        %   theta_d_dot = -k_grad * dJ/d(theta) * sigma
        if abs(w/u) > 0.176 % theta > 10 deg
            % Large dive speeds
            gradJ = -(u * cos(theta) + w * sin(theta));
        else
            % Simplified gradient when theta < 10 deg
            gradJ = -sign(u);
        end
        theta_d = theta_d - h * k_grad * gradJ * sigma;

        % LP filtering the depth command
        z_d = exp(-h*wn_d_z) * z_d + (1 - exp(-h*wn_d_z)) * z_ref;

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

    % MIMO PID controller for pitch and yaw moments
    tau1 = 20; % Constant surge force
    tau5 = -Kp(1,1) * ssa(theta - theta_d) -Kd(1,1) * q - Ki(1,1) * theta_int;
    tau6 = -Kp(2,2) * ssa(psi - psi_d) -Kd(2,2) * r - Ki(2,2) * psi_int;

    tau = [tau1 0 0 0 tau5 tau6]';

    % Depth-attenuated 6-DOF wave forces (compute RAO every 0.1 second)
    if t(i) >= nextRAOtime
        U = sqrt(u^2 + v^2);
        [tau_wave1, waveElevation] = waveForceRAO( ...
            t(i), waveModel, U, psi, environment.beta_wave, zn);
        nextRAOtime = nextRAOtime + RAO_update_period;
    end

    % Ocean current random walks with saturation
    betaVc = sat(betaVc + 0.01 * h * randn, deg2rad(100));
    Vc = sat(Vc + 0.05 * h * randn, 1.0);
    wc = sat(wc + 0.01 * h * randn, 0.2);

    % Store simulation data in a table
    simdata(i,:) = [z_d theta_d psi_d r_d tau' x' Vc betaVc wc ...
        tau_wave1' waveElevation];

    % RK4 method (k+1)
    tau_ext = tau + tau_wave1;
    x = rk4(@lauvMarie, h, x, tau_ext, Vc, betaVc, wc, w0_spectrum);  % LAUV dynamics

    % Euler's integration method (k+1)
    z_int = z_int + h * (zn - z_d);
    theta_int = theta_int + h * ssa(theta - theta_d);
    psi_int = psi_int + h * ssa(psi - psi_d);

end

% ==============================================================================
%% PLOTS
% ==============================================================================
legendLocation = 'best'; legendSize = 11;
if isoctave; legendLocation = 'northeast'; end

% simdata(i,:) = [z_d theta_d psi_d r_d tau' x' Vc betaVc wc ...
%                 tau_wave1' waveElevation];
z_d      = simdata(:,1);
theta_d  = simdata(:,2);
psi_d    = simdata(:,3);
r_d      = simdata(:,4);
tau      = simdata(:,5:10);
nu       = simdata(:,11:16);
eta      = simdata(:,17:22);
Vc       = simdata(:,23);
betaVc   = simdata(:,24);
wc       = simdata(:,25);
tau_wave1 = simdata(:,26:31);
waveElevation = simdata(:,32);
u = nu(:,1); v = nu(:,2); w = nu(:,3);
p = nu(:,4); q = nu(:,5); r = nu(:,6);
phi = eta(:,4); theta = eta(:,5);


% ALOSdata = [y_e z_e alpha_c_hat beta_c_hat]
y_e = ALOSdata(:,1);
z_e = ALOSdata(:,2);
alpha_c_hat = ALOSdata(:,3);
beta_c_hat = ALOSdata(:,4);

uc = Vc .* cos(betaVc);
vc = Vc .* sin(betaVc);

% Crab angles, AOA and SSA
U = sqrt(u.^2+v.^2+w.^2); % Speed
gammaArgument = (u.*sin(theta) - ...
    (v.*sin(phi)+w.*cos(phi)).*cos(theta)) ./ max(U, eps);
gamma = asin(max(-1, min(1, gammaArgument)));
alpha_c = theta - gamma; % Horizontal crab angle
beta_c = atan2(v.*cos(phi)-w.*sin(phi), ...
    u.*cos(theta)+(v.*sin(phi)+w.*cos(phi)).*sin(theta)); % Vertical crab angle
alpha = atan2(w-wc, u-uc); % AOA
beta  = atan2( (v-vc), (u-uc) ); % SSA

%% Generalized velocity
figure(1);
if ~isoctave
    set(gcf, 'Units', 'normalized', ...
        'OuterPosition', [0.01 0.06 0.31 0.88]);
end
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
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Heave position and Euler angles
figure(2);
if ~isoctave
    set(gcf, 'Units', 'normalized', ...
        'OuterPosition', [0.34 0.06 0.31 0.88]);
end
if ControlFlag == 2; z_d = eta(:,3); end
subplot(411),plot(t,eta(:,3),t,z_d,'-.')
xlabel('Time (s)'),title('Down position (m)'),grid
legend('True','Desired')
subplot(412),plot(t,rad2deg(eta(:,4)))
xlabel('Time (s)'),title('Roll angle (deg)'),grid
subplot(413),plot(t,rad2deg(eta(:,5)),t,rad2deg(theta_d),'-.')
xlabel('Time (s)'),title('Pitch angle (deg)'),grid
legend('True','Desired')
subplot(414),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d),'-.')
xlabel('Time (s)'),title('Yaw angle (deg)'),grid
legend('True','Desired')
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Control forces and moments
figure(3);
if ~isoctave
    set(gcf, 'Units', 'normalized', ...
        'OuterPosition', [0.67 0.53 0.31 0.41]);
end
subplot(311),plot(t,tau(:,1))
xlabel('Time (s)'),title('Surge force (N)'),grid
subplot(312),plot(t,tau(:,5))
xlabel('Time (s)'),title('Pitch moment (N m)'),grid
subplot(313),plot(t,tau(:,6))
xlabel('Time (s)'),title('Yaw moment (N m)'),grid
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)

%% Ocean currents and speed
figure(4);
if ~isoctave
    set(gcf, 'Units', 'normalized', ...
        'OuterPosition', [0.67 0.06 0.31 0.41]);
end
subplot(311),plot(t,sqrt(nu(:,1).^2+nu(:,2).^2),t,Vc)
xlabel('Time (s)'),grid
legend('Vehicle horizontal speed (m/s)','Ocean current horizontal speed (m/s)',...
    'Location',legendLocation)
subplot(312),plot(t,nu(:,3),t,wc)
xlabel('Time (s)'),grid
legend('Vehicle heave velocity (m/s)','Ocean current heave velocity (m/s)',...
    'Location',legendLocation)
subplot(313),plot(t,rad2deg(betaVc),'r')
xlabel('Time (s)'),grid
legend('Ocean current direction (deg)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% Crab angles, SSA and AOA
if ControlFlag == 2
    figure(5);
    if ~isoctave
        set(gcf, 'Units', 'normalized', ...
            'OuterPosition', [0.01 0.06 0.31 0.88]);
    end
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
    set(findall(gcf,'type','line'),'linewidth',1.5)
    set(findall(gcf,'type','text'),'FontSize',12)
    set(findall(gcf,'type','legend'),'FontSize',legendSize)
end

%% 2-D position plots with waypoints
if ControlFlag == 2
    figure(6);
    if ~isoctave
        set(gcf, 'Units', 'normalized', ...
            'OuterPosition', [0.33 0.23 0.34 0.50]);
    end
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
    set(findall(gcf, 'type', 'line'), 'LineWidth', 1.5);
    set(findall(gcf, 'type', 'text'), 'FontSize', 12);
    set(findall(gcf, 'type', 'legend'), 'FontSize', legendSize);
end

%% 3-D position plot with waypoints
if ControlFlag == 2
    figure(7);
    if ~isoctave
        set(gcf, 'Units', 'normalized', ...
            'OuterPosition', [0.33 0.23 0.34 0.50]);
    end
    plot3(eta(:,2),eta(:,1),eta(:,3))
    hold on;
    plot3(wpt.pos.y, wpt.pos.x, wpt.pos.z, 'ro', 'MarkerSize', 15);
    hold off
    title('North-East-Down plot (m)')
    xlabel('East'); ylabel('North'); zlabel('Down');
    legend('Actual path','Waypoints','Location',legendLocation),grid
    set(gca, 'ZDir', 'reverse');
    set(findall(gcf,'type','line'),'linewidth',1.5)
    set(findall(gcf,'type','text'),'FontSize',12)
    set(findall(gcf,'type','legend'),'FontSize',legendSize)
    view(-25, 30);  % view(AZ,EL)
end

%% Depth-attenuated 6-DOF first-order wave forces
figure(8); clf;
if ~isoctave
    set(gcf, 'Resize', 'on', 'Units', 'normalized', ...
        'OuterPosition', [0.01 0.06 0.465 0.88]);
end
DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)', ...
    'Roll (N m)', 'Pitch (N m)', 'Yaw (N m)'};
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, tau_wave1(:,DOF), 'LineWidth', 1.5);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF}, 'Location', legendLocation);
end
if ~isoctave
    sgtitle(['Depth-attenuated generalized 1st-order wave forces for ' ...
        '\beta_{wave} = ' num2str(rad2deg(environment.beta_wave)) ...
        '°, H_s = ' num2str(environment.Hs) ' m and \omega_0 = ' ...
        num2str(environment.w0) ' rad/s'], 'FontSize', 12);
end
set(findall(gcf, 'type', 'text'), 'FontSize', 12);
set(findall(gcf, 'type', 'legend'), 'FontSize', legendSize);

%% Directional wave spectrum and wave elevation
figure(9); clf;
subplot(211);
hold on;
if environment.spreadingFlag
    plot(environment.Omega, ...
        environment.S_M(:,floor(length(environment.mu)/2)), ...
        'b->', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(environment.Omega, ...
        environment.S_M(:,floor(length(environment.mu)/4)), ...
        'k-o', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(environment.Omega, environment.S_M(:,length(environment.mu)), ...
        'g-', 'LineWidth', 2.0);
    plot([environment.w0 environment.w0], ...
        [min(environment.S_M(:)) max(environment.S_M(:))], ...
        'r-.', 'LineWidth', 1.5);
    legend('\mu = 0 deg', '\mu = 45 deg', '\mu = 90 deg', ...
        ['\omega_0 = ' num2str(environment.w0) ' rad/s'], ...
        'Location', legendLocation);
else
    plot(environment.Omega, environment.S_M(:,1), ...
        'b-', 'LineWidth', 1.5);
    plot([environment.w0 environment.w0], ...
        [min(environment.S_M(:)) max(environment.S_M(:))], ...
        'r-.', 'LineWidth', 1.5);
    legend('S(\omega)', ...
        ['\omega_0 = ' num2str(environment.w0) ' rad/s'], ...
        'Location', legendLocation);
end
hold off;
xlabel('Omega (rad/s)');
ylabel('m^2 s');
title([environment.spectrumType ' spectrum']);
grid on;

subplot(212);
plot(t, waveElevation, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('m');
title(['Wave elevation for wave direction \beta_{wave} = ' ...
    num2str(rad2deg(environment.beta_wave)) '°, H_s = ' ...
    num2str(environment.Hs) ' m and \omega_0 = ' ...
    num2str(environment.w0) ' rad/s']);
grid on;
set(findall(gcf, 'type', 'text'), 'FontSize', 12);
set(findall(gcf, 'type', 'legend'), 'FontSize', legendSize);

% Display the vehicle data and a photograph of LAUV Marie at NTNU.
vehicleData = {...
    'Length', sprintf('%.2f m', vessel.main.Lpp), ...
    'Diameter', sprintf('%.2f m', vessel.main.B), ...
    'Mass', sprintf('%.1f kg', vessel.main.m), ...
    'Reference depth', sprintf('%.1f m', vessel.main.submergenceDepth), ...
    'Hydrodynamic model', 'MSS-Capytaine', ...
    'Photo credit', 'NTNU AUR-Lab.'};
displayVehicleData('Light AUV (LAUV) Marie', vehicleData, ...
    'lauv_marie.png', 10);


end

%% FUNCTIONS
function [ControlFlag, Hs, w0_spectrum, beta_wave, ...
    Vc, betaVc, wc] = controlMethod()

f = figure('Position', [300, 100, 900, 650], ...
    'Name', 'Simulation Parameters', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'Resize', 'off', ...
    'WindowStyle', 'modal');

% Add button group for control methods
bg1 = uibuttongroup('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [15 475 870 160], ...
    'Title', 'Control Methods', ...
    'FontSize',14, ...
    'FontWeight','bold');
radio1 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'PID pole-placement for heading and gradient-based depth control', ...
    'Position', [15 75 825 35], ...
    'Tag', '1', ...
    'Value', 1);
radio2 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize',13, 'String', ...
    'ALOS guidance law for 3-D path following', ...
    'Position', [15 30 825 35], ...
    'Tag', '2');

% Wave-spectrum parameters
bg2 = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [15 285 870 165], ...
    'Title', 'Directional JONSWAP Wave Spectrum', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

uicontrol(bg2, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Significant wave height, Hs (m)', ...
    'Position', [20 85 250 30]);
editHs = uicontrol(bg2, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '0.3', ...
    'Position', [20 45 200 34]);

uicontrol(bg2, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Peak frequency, w0 (rad/s)', ...
    'Position', [300 85 250 30]);
editW0 = uicontrol(bg2, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '2.7', ...
    'Position', [300 45 200 34]);

uicontrol(bg2, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Mean wave direction (deg, NED)', ...
    'Position', [580 85 250 30]);
editBetaWave = uicontrol(bg2, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '140', ...
    'Position', [580 45 200 34]);

% Ocean-current parameters
bg3 = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [15 95 870 165], ...
    'Title', 'Ocean Currents', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

uicontrol(bg3, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Horizontal speed, Vc (m/s)', ...
    'Position', [20 85 250 30]);
editVc = uicontrol(bg3, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '0.3', ...
    'Position', [20 45 200 34]);

uicontrol(bg3, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Direction, betaVc (deg NED)', ...
    'Position', [300 85 280 30]);
editBetaVc = uicontrol(bg3, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '20', ...
    'Position', [300 45 200 34]);

uicontrol(bg3, ...
    'Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 13, ...
    'String', 'Vertical speed, wc (m/s, +down)', ...
    'Position', [580 85 280 30]);
editWc = uicontrol(bg3, ...
    'Style', 'edit', ...
    'FontSize', 13, ...
    'String', '0.1', ...
    'Position', [580 45 200 34]);

% Add OK button to confirm selections
uicontrol('Style', 'pushbutton', ...
    'String', 'OK', ...
    'FontSize', 13, ...
    'Position', [15 30 110 42], ...
    'Callback', @(src, evt) uiresume(f));

uiwait(f); % Wait for uiresume to be called on figure handle

% Determine which control method was selected
if get(radio1, 'Value') == 1
    ControlFlag = str2double(get(radio1, 'Tag'));
else
    ControlFlag = str2double(get(radio2, 'Tag'));
end

% Read and validate wave-spectrum and ocean-current parameters
Hs = str2double(get(editHs, 'String'));
w0_spectrum = str2double(get(editW0, 'String'));
betaWaveDeg = str2double(get(editBetaWave, 'String'));
Vc = str2double(get(editVc, 'String'));
betaVcDeg = str2double(get(editBetaVc, 'String'));
wc = str2double(get(editWc, 'String'));

if ~isscalar(Hs) || ~isfinite(Hs) || Hs <= 0
    close(f);
    error('Significant wave height Hs must be a positive finite scalar.');
end
if ~isscalar(w0_spectrum) || ~isfinite(w0_spectrum) || w0_spectrum <= 0
    close(f);
    error('Peak frequency w0 must be a positive finite scalar.');
end
if ~isscalar(betaWaveDeg) || ~isfinite(betaWaveDeg)
    close(f);
    error('Mean wave direction must be a finite scalar in degrees.');
end
if ~isscalar(Vc) || ~isfinite(Vc) || Vc < 0
    close(f);
    error('Horizontal current speed Vc must be finite and nonnegative.');
end
if ~isscalar(betaVcDeg) || ~isfinite(betaVcDeg)
    close(f);
    error('Horizontal current direction betaVc must be finite in degrees.');
end
if ~isscalar(wc) || ~isfinite(wc)
    close(f);
    error('Vertical current speed wc must be a finite scalar.');
end
beta_wave = deg2rad(betaWaveDeg);
betaVc = deg2rad(betaVcDeg);

close(f);  % close the figure after obtaining the selections

disp('-----------------------------------------------------------------------');
disp('MSS toolbox: LAUV Marie AUV');
disp('Kinematics:          Euler angle representation (12 states)');

if (ControlFlag == 1)
    disp('Heading autopilot:   MIMO PID pole-placement control');
    disp('Depth autopilot:     Gradient-based depth control');
else
    disp('Attitude autopilot:  MIMO PID pole-placement control');
    disp('Path-following:      ALOS guidance law for 3-D path following');
end
disp('------------------------------------------------------------------------');
disp('Simulating...');

end
