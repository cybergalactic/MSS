function SIMotter()
% SIMOTTER Unified 6-DOF maneuvering and seakeeping simulation of the Otter USV.
% The default otter.m model combines MSS-Capytaine power-based equivalent
% added mass, radiation damping and hydrostatic restoring with rigid-body
% and payload dynamics, nonlinear damping, twin-pontoon cross-flow drag,
% Coriolis and centripetal terms, propeller thrust, and ocean currents.
% First-order wave forces and moments are synthesized from the Capytaine
% force RAOs and a directional JONSWAP spectrum. The loads are transformed
% from the dataset CG to CO and applied through the selected mass matrix at
% every RK4 stage, with RAO updates every 0.1 s. The RAOs use the nominal 80 kg
% load condition for either dynamics model. waveForceRAO accumulates signed
% encounter phase during changes in speed and heading, as in Fossen (2027),
% Chapter 10. Wave-force histories, the spectrum and wave elevation are plotted.
% Default sea-state parameters are Hs = 0.20 m and w0 = 3.0 rad/s (peak period 2.09 s)
%
% Guidance and control options include heading control and path following
% along straight-line segments or cubic Hermite splines.
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
%   otter.m                 - MSS-Capytaine Otter dynamics (default)
%   otter_legacy.m          - Legacy first-principles Otter model
%   refModel.m              - Reference model for autopilot systems
%   ILOSpsi.m               - ILOS guidance for waypoint path following
%   ALOSpsi.m               - ALOS guidance for waypoint path following
%   ALOSpsiHermite          - ALOS guidance using cubic Hermite splines
%   hermiteSpline.m         - Cubic Hermite spline computations
%   waveInitialization.m    - Directional spectrum and force-RAO setup
%   waveForceRAO.m          - First-order wave forces and elevation
%   LOSobserver.m           - Observer for LOS guidance
%
% The setup GUI and configuration validation are local functions below.
%
% References:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%     Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
%   T. I. Fossen and A. P. Aguiar (2024). A Uniform Semiglobal Exponential
%     Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path Following.
%     Automatica, 163, 111556. https://doi.org/10.1016/j.automatica.2024.111556
%
%   T. I. Fossen (2023). An Adaptive Line-of-sight (ALOS) Guidance Law for
%     Path Following of Aircraft and Marine Craft. IEEE Transactions on Control
%     Systems Technology, 31(6), 2887-2894.
%     https://doi.org/10.1109/TCST.2023.3259819
%
% Author: Thor I. Fossen
% Date: 2021-04-25
% Revisions:
%   2023-10-14: Added heading autopilot and reference model.
%   2024-04-01: Integrated ALOS/ILOS path-following control algorithms for
%               straight-line paths and Hermite splines.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method with RK4.
%   2026-06-27: Added propeller speed saturation using limits from otter.m.
%   2026-10-07: Corrected reverse-thrust control allocation (E. Krizman)
%   2026-10-18: Smoothed Hermite turns using shared waypoint tangents and
%               chord-length spacing; kept intermediate waypoints only for
%               straight-line guidance and corrected Hermite search indexing.
%   2026-10-26: Added an optional legacy-model checkbox to the control GUI.
%               Default to the MSS-Capytaine model; use otterDynamics for
%               both controller initialization and RK4 propagation.
%               Added current, wave and simulation settings to the setup GUI.
%               Accumulate encounter phase in waveForceRAO as in Chapter 10.

clearvars;                                  % Clear variables from memory
close all;                                  % Close all figure windows
clear ALOSpsi ILOSpsi ALOSpsiHermite        % Clear persistent variables

%% SIMULATION SETUP
cfg = guiSIMotter();
if isempty(cfg)
    return;
end
h = cfg.simulation.h;
T_final = cfg.simulation.T_final;
environment = cfg.environment;
V_c = environment.Vc;
beta_c = environment.betaVc;
ControlFlag = cfg.controlFlag;
otterDynamics = @otter;
if cfg.useLegacy
    otterDynamics = @otter_legacy;
end

% First-order wave excitation: the Capytaine RAOs are referenced to the
% dataset CG. Transform their generalized forces to the dynamics origin CO.
RAO_update_period = 0.1;
nextRAOtime = 0;
tau_wave1 = zeros(6,1);
waveElevation = 0;
waveModel = [];
H_wave = eye(6);
if environment.wavesEnabled && environment.Hs > 0
    data = load(resolveHydroVesselFile('otter_80kg.mat'), 'vessel');
    H_wave = Hmtrx(data.vessel.main.CG(:));
    rng(1); % Reproducible stochastic wave realization
    [environment, waveModel] = waveInitialization(data.vessel, environment);
end

% Load condition
mp = 25;                   % Payload mass (kg), maximum value 45 kg
r_bP = [0.19 0 -0.35]';    % Payload position (m), about -1.5 deg pitch at 100 N

% Waypoints
wpt.pos.x = [0 0   150 150 -100 -100 200];
wpt.pos.y = [0 200 200 -50  -50  250 250];

% ALOS and ILOS parameters
Delta_h = 10;                    % Look-ahead distance
gamma_h = 0.001;                 % ALOS adaptive gain
kappa_int = 0.001;               % ILOS integral gain

% Additional parameter for straight-line path following
R_switch = 5;                    % Radius of switching circle
K_f = 0.3;                       % LOS observer gain

% Initial heading, vehicle points towards next waypoint
psi0 = atan2(wpt.pos.y(2) - wpt.pos.y(1), wpt.pos.x(2) - wpt.pos.x(1));

% Additional parameters for Hermite spline path following
Umax = 3;                        % Maximum speed for Hermite spline LOS
idx_start = 1;                   % Initial index for Hermite spline
[w_path, x_path, y_path, dx, dy, pi_h, pp_x, pp_y, N_horizon, ~] = ...
    hermiteSpline(wpt, Umax, h); % Compute Hermite spline for path following

% Otter USV input matrix
[~,~,M,B_prop,n_min,n_max] = otterDynamics();
% The /2 converts the combined twin-propeller coefficients to per-propeller
% values. These coefficients should match the selected Otter model.
k_pos = 0.02216 / 2;
k_neg = 0.01289 / 2;
% B_prop maps n.*abs(n) through k_pos. Dividing by k_pos gives the
% geometry-only matrix mapping individual propeller thrusts to [tau_X; tau_N].
B_thrust = B_prop / k_pos;
B_thrust_inv = invQR(B_thrust);

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

% Display simulation options
if ControlFlag ~= 4
    % Preserve intermediate waypoints for straight-line guidance only.
    wpt = addIntermediateWaypoints(wpt, 2);
end
displayControlMethod(ControlFlag, R_switch, Delta_h, environment);

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
simdata = zeros(nTimeSteps, 21);    % Preallocate table for simulation data

for i = 1:nTimeSteps

    % Measurements
    r = x(6);       % Yaw rate
    xn = x(7);      % North position
    yn = x(8);      % East position
    psi = x(12);    % Yaw angle

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
            psi_ref = ILOSpsi(xn, yn, Delta_h, kappa_int, h, R_switch, wpt);
            [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

        case 4  % ALOS heading autopilot, cubic Hermite spline interpolation
            [psi_ref, idx_start] = ALOSpsiHermite(w_path, x_path, ...
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
    thrust_c = B_thrust_inv * [tau_X; tau_N]; % Required propeller thrusts
    n_c = zeros(2,1);
    for j = 1:2
        if thrust_c(j) >= 0
            n_c(j) = sqrt(thrust_c(j) / k_pos);
        else
            n_c(j) = -sqrt(-thrust_c(j) / k_neg);
        end
    end
    n_c = satlim(n_c, n_min, n_max);  % Saturate propeller command

    % First-order surface-wave loads, held between RAO updates
    if ~isempty(waveModel) && t(i) + 10*eps(max(1,t(i))) >= nextRAOtime
        U = hypot(x(1), x(2));
        [tau_CG, waveElevation, waveModel] = waveForceRAO(t(i), waveModel, ...
            U, x(12), environment.beta_wave);
        tau_wave1 = H_wave' * tau_CG;

        nextRAOtime = (floor((t(i) + 10*eps(max(1,t(i)))) / ...
            RAO_update_period) + 1) * RAO_update_period;
    end

    % Store states, references, wave forces at CO, and wave elevation.
    simdata(i, :) = [x', r_d, psi_d, tau_wave1', waveElevation];

    % RK4 method x(k+1)
    x = rk4(@otterWithWaves, h, x, otterDynamics, n, mp, r_bP, ...
        V_c, beta_c, tau_wave1);

    % Propeller dynamics n(k+1)
    n = n + h/T_n * (n_c - n);     % Euler's method
    n = satlim(n, n_min, n_max);   % Saturate actual propeller speed after update

    % Heading integral state z(k+1)
    z_psi = z_psi + h * ssa(psi - psi_d); % Euler's method

end

% ------------------------------------------------------------------------------
%% PLOTS
% ------------------------------------------------------------------------------
scrSz = get(0, 'ScreenSize'); % Get screen dimensions
legendSize = 12;

% Simulation data structure
nu   = simdata(:,1:6);
eta  = simdata(:,7:12);
r_d = simdata(:,13);
psi_d = simdata(:,14);
tau_wave1 = simdata(:,15:20);
waveElevation = simdata(:,21);

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

xlabel('East (m)', 'FontSize', 12);
ylabel('North (m)', 'FontSize', 12);
title('North-East Positions (m)', 'FontSize', 12);
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
set(findall(gcf,'type','text'),'FontSize',12);
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
set(findall(gcf,'type','text'),'FontSize',12);
set(findall(gcf,'type','legend'),'FontSize',legendSize);

% First-order wave forces and spectrum, as in SIMlauvMarie.
if ~isempty(waveModel)
    legendLocation = 'best';
    if isoctave; legendLocation = 'northeast'; end
    figure(5); clf;
    DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)', ...
        'Roll (N m)', 'Pitch (N m)', 'Yaw (N m)'};
    for DOF = 1:6
        subplot(6,1,DOF);
        plot(t, tau_wave1(:,DOF), 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel(DOF_txt{DOF}); grid on;
    end
    if ~isoctave
        sgtitle('Generalized first-order wave forces at CO');
    end
    set(findall(gcf,'type','text'),'FontSize',12);

    figure(6); clf;
    subplot(211); hold on;
    if environment.spreadingFlag
        % Use the actual randomized spreading angles in the legend.
        angles = [0 pi/4 pi/2];
        for k = 1:length(angles)
            [~, idx] = min(abs(environment.mu - angles(k)));
            plot(environment.Omega, environment.S_M(:,idx), ...
                'LineWidth', 1.5, 'DisplayName', ...
                sprintf('\\mu = %.1f deg', rad2deg(environment.mu(idx))));
        end
        ylabel('m^2 s/rad');
    else
        plot(environment.Omega, environment.S_M(:,1), ...
            'LineWidth', 1.5, 'DisplayName', 'S(\omega)');
        ylabel('m^2 s');
    end
    plot([environment.w0 environment.w0], ...
        [0 max(environment.S_M(:))], 'r-.', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\omega_0 = %.2f rad/s', environment.w0));
    hold off; xlabel('Omega (rad/s)');
    title([environment.spectrumType ' spectrum']); grid on;
    legend('show', 'Location', legendLocation);
    subplot(212);
    plot(t, waveElevation, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('m'); grid on;
    title(sprintf('Wave elevation: H_s = %.2f m, \\beta_{wave} = %.1f deg', ...
        environment.Hs, rad2deg(environment.beta_wave)));
    set(findall(gcf,'type','text'),'FontSize',12);
end

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

% ------------------------------------------------------------------------------
%% HELPER FUNCTION FOR PLOTTING
% ------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
%% DISPLAY CONTROL METHOD
% ------------------------------------------------------------------------------
function displayControlMethod(ControlFlag, R_switch, Delta_h, environment)
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
            disp('ALOS path-following control using cubic Hermite splines');
            disp(['Look-ahead distance: Delta_h = ', num2str(Delta_h), ' m']);
    end
    fprintf('Ocean current: %.2f m/s, direction %.1f deg (towards)\n', ...
        environment.Vc, rad2deg(environment.betaVc));
    if environment.wavesEnabled
        fprintf(['Sea state: Hs = %.2f m, peak frequency = %.2f rad/s, ' ...
            'direction = %.1f deg (towards)\n'], ...
            environment.Hs, environment.w0, rad2deg(environment.beta_wave));
        disp('Wave excitation: Capytaine first-order force RAOs (both models).');
    else
        disp('Sea state: calm water (waves off)');
    end
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end


% ------------------------------------------------------------------------------
%% LOCAL SETUP GUI
% ------------------------------------------------------------------------------
function cfg = guiSIMotter()
% GUISIMOTTER Select dynamics, control, current, waves and simulation timing.
% cfg = guiSIMotter() returns validated SIMotter settings. Closing the
% window or choosing Cancel returns []. GUI directions are in degrees;
% cfg.environment stores directions in radians, towards in NED axes.
%
% Author: Thor I. Fossen
% Date:   2026-10-26

defaults = otterSimulationConfig();
cfg = [];
methods = {'PID heading autopilot, no path following', ...
    'ALOS path following: straight lines and waypoint switching', ...
    'ILOS path following: straight lines and waypoint switching', ...
    'ALOS path following: cubic Hermite splines'};

f = figure('Name', 'Otter Simulation Parameters', ...
    'Position', [300 80 900 780], 'MenuBar', 'none', ...
    'NumberTitle', 'off', 'Resize', 'off', 'WindowStyle', 'modal', ...
    'CloseRequestFcn', @onCancel);

bgControl = uibuttongroup('Parent', f, 'Units', 'pixels', ...
    'Position', [15 555 870 210], 'Title', 'Control Methods', ...
    'FontSize', 14, 'FontWeight', 'bold');
controlRadio = cell(1, length(methods));
for k = 1:length(methods)
    controlRadio{k} = uicontrol(bgControl, 'Style', 'radiobutton', ...
        'String', methods{k}, 'FontSize', 13, 'Tag', num2str(k), ...
        'Position', [20 145-(k-1)*45 825 35]);
end
set(bgControl, 'SelectedObject', controlRadio{defaults.controlFlag});

pWave = uipanel('Parent', f, 'Units', 'pixels', ...
    'Position', [15 350 870 190], 'Title', 'Directional JONSWAP Wave Spectrum', ...
    'FontSize', 14, 'FontWeight', 'bold');
wavesCheckbox = uicontrol(pWave, 'Style', 'checkbox', ...
    'String', 'Enable waves', 'FontSize', 13, ...
    'Value', defaults.environment.wavesEnabled, 'Position', [20 125 825 30]);
editHs = panelEdit(pWave, 'Significant wave height, Hs (m)', 20, 75, ...
    defaults.environment.Hs);
editW0 = panelEdit(pWave, 'Peak frequency, w0 (rad/s)', 300, 75, ...
    defaults.environment.w0);
editBetaWave = panelEdit(pWave, 'Wave direction, towards (deg)', 580, 75, ...
    rad2deg(defaults.environment.beta_wave));

pCurrent = uipanel('Parent', f, 'Units', 'pixels', ...
    'Position', [15 190 870 145], 'Title', 'Ocean Currents and Dynamics', ...
    'FontSize', 14, 'FontWeight', 'bold');
editVc = panelEdit(pCurrent, 'Horizontal speed, Vc (m/s)', 20, 75, ...
    defaults.environment.Vc);
editBetaVc = panelEdit(pCurrent, 'Current direction, towards (deg)', 300, 75, ...
    rad2deg(defaults.environment.betaVc));
legacyCheckbox = uicontrol(pCurrent, 'Style', 'checkbox', ...
    'String', 'Legacy Otter (2019-2026)', 'FontSize', 12, ...
    'Value', defaults.useLegacy, 'Position', [580 75 270 30]);
addText(pCurrent, {'Unchecked: unified MSS-Capytaine', ...
    'model with nonlinear damping.'}, [580 20 270 48], 11);
addText(pCurrent, 'Directions: towards, clockwise from North; East = 90 deg.', ...
    [20 10 540 23], 11);

pSimulation = uipanel('Parent', f, 'Units', 'pixels', ...
    'Position', [15 65 870 110], 'Title', 'Simulation', ...
    'FontSize', 14, 'FontWeight', 'bold');
editDuration = panelEdit(pSimulation, 'Duration, T_final (s)', 20, 50, ...
    defaults.simulation.T_final);
editStep = panelEdit(pSimulation, 'Time step, h (s)', 300, 50, ...
    defaults.simulation.h);
addText(pSimulation, {'Time step: 0 < h <= 0.1 s.', ...
    'Default: h = 0.05 s.'}, [580 15 270 60], 12);

uicontrol(f, 'Style', 'pushbutton', 'String', 'Start simulation', ...
    'Position', [15 12 170 40], 'FontSize', 13, 'Callback', @onStart);
uicontrol(f, 'Style', 'pushbutton', 'String', 'Cancel', ...
    'Position', [200 12 110 40], 'FontSize', 13, 'Callback', @onCancel);

drawnow;
uiwait(f);
if ishandle(f)
    delete(f);
end

    function onStart(~, ~)
        selected = defaults;
        selected.controlFlag = str2double(get(get(bgControl, 'SelectedObject'), 'Tag'));
        selected.useLegacy = get(legacyCheckbox, 'Value');
        selected.simulation.T_final = str2double(get(editDuration, 'String'));
        selected.simulation.h = str2double(get(editStep, 'String'));
        selected.environment.Vc = str2double(get(editVc, 'String'));
        selected.environment.betaVc = deg2rad(str2double(get(editBetaVc, 'String')));
        selected.environment.wavesEnabled = get(wavesCheckbox, 'Value');
        selected.environment.Hs = str2double(get(editHs, 'String'));
        selected.environment.w0 = str2double(get(editW0, 'String'));
        selected.environment.beta_wave = deg2rad(str2double(get(editBetaWave, 'String')));
        try
            selected = otterSimulationConfig(selected);
        catch err
            errordlg(err.message, 'Check simulation settings');
            return;
        end
        cfg = selected;
        uiresume(f);
    end

    function onCancel(~, ~)
        cfg = [];
        uiresume(f);
    end
end

function handle = panelEdit(parent, label, x, y, value)
addText(parent, label, [x y 275 30], 13);
handle = uicontrol(parent, 'Style', 'edit', 'String', num2str(value,12), ...
    'FontSize', 13, 'Position', [x y-40 200 34]);
end

function addText(parent, text, position, fontSize)
if nargin < 4
    fontSize = 10;
end
uicontrol('Parent', parent, 'Style', 'text', 'String', text, ...
    'HorizontalAlignment', 'left', 'FontSize', fontSize, 'Position', position);
end


% ------------------------------------------------------------------------------
%% LOCAL CONFIGURATION AND VALIDATION
% ------------------------------------------------------------------------------
function cfg = otterSimulationConfig(cfg)
% OTTERSIMULATIONCONFIG Default settings and validation for SIMotter.
% cfg = otterSimulationConfig() returns the default simulation environment.
% cfg = otterSimulationConfig(cfg) validates settings and wraps directions
% to [0,2*pi). Directions are towards, clockwise from North, in radians.
%
% Author: Thor I. Fossen
% Date:   2026-10-26

if nargin == 0
    cfg.controlFlag = 1;
    cfg.useLegacy = false;
    cfg.simulation.T_final = 700;
    cfg.simulation.h = 0.05;
    cfg.environment.Vc = 0.3;
    cfg.environment.betaVc = deg2rad(30);
    cfg.environment.wavesEnabled = false;
    cfg.environment.Hs = 0.2;
    cfg.environment.w0 = 3.0;
    cfg.environment.beta_wave = deg2rad(140);
    % Directional spectrum used by waveInitialization / waveForceRAO.
    cfg.environment.spectrumType = 'JONSWAP';
    cfg.environment.spreadingFlag = true;
    cfg.environment.numFreqIntervals = 100;
    cfg.environment.numDirections = 24;
    return;
end

checkScalar(cfg.controlFlag, 'Control method');
if cfg.controlFlag ~= fix(cfg.controlFlag) || ...
        cfg.controlFlag < 1 || cfg.controlFlag > 4
    error('SIMotter:ControlMethod', 'Select one of the four control methods.');
end
checkSwitch(cfg.useLegacy, 'Legacy-model selection');
checkScalar(cfg.simulation.T_final, 'Simulation duration');
checkScalar(cfg.simulation.h, 'Time step');
if cfg.simulation.T_final <= 0
    error('SIMotter:Duration', 'Simulation duration must be positive.');
end
if cfg.simulation.h <= 0 || cfg.simulation.h > 0.1
    error('SIMotter:TimeStep', 'Time step must be positive and no greater than 0.1 s.');
end
if cfg.simulation.h > cfg.simulation.T_final
    error('SIMotter:Duration', 'Simulation duration must be at least one time step.');
end

checkScalar(cfg.environment.Vc, 'Current speed');
checkScalar(cfg.environment.betaVc, 'Current direction');
checkSwitch(cfg.environment.wavesEnabled, 'Wave selection');
checkScalar(cfg.environment.Hs, 'Significant wave height');
checkScalar(cfg.environment.w0, 'Peak frequency');
checkScalar(cfg.environment.beta_wave, 'Wave direction');
if cfg.environment.Vc < 0
    error('SIMotter:CurrentSpeed', 'Current speed must be nonnegative.');
end
if cfg.environment.Hs < 0
    error('SIMotter:WaveHeight', 'Significant wave height must be nonnegative.');
end
if cfg.environment.w0 <= 0
    error('SIMotter:PeakFrequency', 'Peak frequency must be positive.');
end
cfg.useLegacy = logical(cfg.useLegacy);
cfg.environment.wavesEnabled = logical(cfg.environment.wavesEnabled);
cfg.environment.betaVc = mod(cfg.environment.betaVc, 2*pi);
cfg.environment.beta_wave = mod(cfg.environment.beta_wave, 2*pi);
end

function checkScalar(value, label)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
    error('SIMotter:Scalar', '%s must be a finite real number.', label);
end
end

function checkSwitch(value, label)
if ~(isnumeric(value) || islogical(value)) || ~isreal(value) || ...
        ~isscalar(value) || ~isfinite(value) || ~any(value == [0 1])
    error('SIMotter:Switch', '%s must be on or off.', label);
end
end

function xdot = otterWithWaves(x, dynamics, n, mp, r_bP, V_c, beta_c, tau_wave)
% Apply external forces at CO through the selected model's mass matrix.
% Evaluate dynamics at every RK4 stage; the wave force is held between updates.
[xdot, ~, M] = dynamics(x, n, mp, r_bP, V_c, beta_c);
xdot(1:6) = xdot(1:6) + M \ tau_wave;
end
