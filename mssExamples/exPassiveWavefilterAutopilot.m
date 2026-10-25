% Passive wave filtering and PID heading control of the Mariner-class cargo ship.
% Select the measured motion signals in the popup window:
%
%   'Yaw angle only': Estimate LF heading, yaw rate and rudder bias.
%          Measured rudder angle is a known input to the Nomoto model.
%   'Yaw angle and yaw rate': Correct gyro bias, then filter
%      the measured yaw rate. No Nomoto model is needed by this observer.
%
% Both observers use e = ssa(y - psi_hat - psi_w_hat) and
%
%   xi_w_hat_dot  = psi_w_hat + k1(1)*e
%   psi_w_hat_dot = -w_0^2*xi_w_hat - 2*lambda*w_0*psi_w_hat + k1(2)*e
%   psi_hat_dot   = r_hat + k2*e
%
% Yaw angle only: r_hat_dot = (-r_hat + K*(delta + b_hat))/T + k3*e,
%            b_hat_dot = -b_hat/T_b + k4*e.
%            Positive rudder bias adds to the rudder input, as in the book.
% Yaw angle and yaw rate: b_hat_dot = -b_hat/T_b + k3*e; b_hat is rate-gyro bias (rad/s).
%            r_hat = LPF(notchFilter(y_r)) + b_hat,
%            LPF(s) = (w_lp/(s + w_lp))^2 (two first-order stages).
% Both passive observers retain the bias decay term required by the design.
% Finite T_b allows a residual bias error for a constant disturbance.
%
% Linear observer dynamics are discretized with expm and held measurements.
% The gyro notch uses notchFilter(...,'RK4'), matching the transfer
% function (s^2 + 2*zeta_n*w_0*s + w_0^2)/(s^2 + 2*w_0*s + w_0^2).
% Uses base MATLAB and MSS headingAutopilot, ssa, sat, notchFilter and lowPassFilter.
% The ship, waves and PID reference model follow exKalmanWavefilterAutopilot.m.
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-10-25

clearvars; close all;
clear headingAutopilot;          % Reset the PID and reference-model states
rng(23341);                      % Repeatable wave realization

% ------------------------------------------------------------------------------
%% USER INPUTS
% ------------------------------------------------------------------------------
sensor_option = selectMeasurements;
if sensor_option == 0
    return;                      % Stop if the selection window is closed
end
b_gyro = deg2rad(0.05);          % Unknown rate-gyro bias, with rate-gyro measurements (rad/s)
h = 0.1;                         % Sampling time (s)
T_final = 600;                   % Simulation time (s)
psi_ref = deg2rad(30.0);         % Commanded heading (rad)
b_true = deg2rad(1.0);           % Unknown equivalent rudder bias (rad)

% Mariner ship and rudder servo
K = 0.185;
T1 = 118; T2 = 7.8; T3 = 18.5;   % Second-order Nomoto time constants (s)
T = T1 + T2 - T3;                % First-order Nomoto time constant (s)
T_delta = 1;                     % Rudder servo time constant (s)
delta_max = deg2rad(20);         % Rudder angle limit (rad)
rate_max = deg2rad(5);           % Rudder rate limit (rad/s)

% MSS PID heading autopilot and its internal reference model
wn = 0.1; zeta = 0.5;            % Td = 10 s and Ti = 100 s, as in the case study
wn_d = 0.05;                     % Reference-model bandwidth (rad/s)
r_max = deg2rad(1);              % Desired yaw-rate limit (rad/s)
M_nomoto = T/K;                  % Nomoto inertia for a rudder-angle input

% Wave model: deg2rad(3) * 2*lambda*w_0*sigma*s / (s^2 + 2*lambda*w_0*s + w_0^2)
w_0 = 1.2; lambda = 0.1; sigma = 0.5;
wave_gain = deg2rad(3);
Aw = [0 1; -w_0^2 -2*lambda*w_0];
Bw = [0; wave_gain * 2*lambda*w_0*sigma];
Mw = expm(h * [Aw Bw; zeros(1,3)]);
Adw = Mw(1:2,1:2); Bdw = Mw(1:2,3);

% Ship states: [psi r r_dot delta]'. Input: held rudder rate delta_dot.
As = [0 1 0 0
      0 0 1 0
      0 -1/(T1*T2) -(T1+T2)/(T1*T2) K/(T1*T2)
      0 0 0 0];
Bs = [0 0 K*T3/(T1*T2) 1]';
Es = [0 0 K/(T1*T2) 0]';       % Book convention: delta + b_true
Ms = expm(h * [As Bs Es; zeros(2,6)]);
Ads = Ms(1:4,1:4); Bds = Ms(1:4,5); Eds = Ms(1:4,6);

% Passive observer gains and LF models (Chapter 13).
% Both cases use [xi_w_hat psi_w_hat psi_hat r_filtered b_hat]'.
k = zeros(5,1);
switch sensor_option
    case 1                     % Yaw angle measured; yaw rate estimated
        sensor_name = 'Yaw angle only';
        T_b = 20000; w_c = 1.3; k3 = 0.04; k4 = 0.3;
        % Passive-gain constraint: 0 < 1/T_b < k3/k4 < w_0 < w_c.
        assert(1/T_b < k3/k4 && k3/k4 < w_0 && w_0 < w_c, ...
            'The passive observer gain constraint must be satisfied.');
        rate_model = [-1/T K/T];
        bias_in_heading = 0;
        B = [0 0 0 K/T 0]';
        k(4:5) = [k3; k4];
    case 2                     % Yaw angle and yaw rate measured
        sensor_name = 'Yaw angle and yaw rate';
        T_b = 10000; w_c = 1.3; k3 = 0.02;
        w_lp = 0.5;            % Bandwidth of each gyro low-pass stage (rad/s)
        rate_model = [-w_lp 0];
        bias_in_heading = 1;
        B = [0 0 0 w_lp 0]';
        k(5) = k3;
end
k1 = [-2*w_c*(1-lambda)/w_0; 2*w_0*(1-lambda)];
k2 = w_c;
k(1:3) = [k1; k2];
C = [0 1 1 0 0];
A_LF = [0 1 bias_in_heading; 0 rate_model; 0 0 -1/T_b];
A = blkdiag(Aw, A_LF);          % Common WF model and sensor-specific LF model
PHI = expm(h * [A-k*C k B; zeros(2,7)]);
Ad = PHI(1:5,1:5); Bd = PHI(1:5,6:7);
zeta_n = 0.01;                  % Notch attenuation at w_0
xf = zeros(2,1);                % notchFilter RK4 states
xf_lp = 0;                     % First gyro low-pass stage
x_hat = zeros(5,1);
xs = zeros(4,1);                % Ship states [psi r r_dot delta]'
xw = zeros(2,1);                % Wave states [xi_w psi_w]'
t = (0:h:T_final)';
nTimeSteps = length(t);
% Columns: psi, r, delta, delta_c, psi_w, y, psi_w_hat, psi_hat, r_hat,
%          estimated bias, true bias, measured yaw rate (with rate-gyro measurements only).
simdata = zeros(nTimeSteps,12);
wave_noise = randn/sqrt(h);

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
for i = 1:nTimeSteps

    % Compass measures LF heading + WF heading
    psi = xs(1); r = xs(2); delta = xs(4);
    psi_w = xw(2);
    y = psi + psi_w;
    e = ssa(y - C*x_hat);

    switch sensor_option
        case 1                 % Reconstruct yaw rate from compass and rudder
            r_hat = x_hat(4); b_hat = x_hat(5);
            bias_true = b_true; y_r = NaN;
            u_obs = delta;     % Measured rudder angle is a known input
        case 2                 % Filter the gyro measurement; correct its bias
            r_w = Aw(2,:)*xw + Bw(2)*wave_noise;
            y_r = r + r_w - b_gyro; % Book convention: u_ars = r - b_ars
            b_hat = x_hat(5); r_hat = x_hat(4) + b_hat;
            bias_true = b_gyro;
            [xf, u_obs] = notchFilter(xf, y_r, w_0, zeta_n, h, 'RK4');
            [xf_lp, u_obs] = lowPassFilter(xf_lp, u_obs, w_lp, h);
    end

    delta_c = headingAutopilot(x_hat(3),r_hat,psi_ref, ...
        r_max,M_nomoto,wn,zeta,h,wn_d);
    delta_ref = sat(delta_c, delta_max);
    delta_dot = sat((delta_ref-delta)/T_delta, rate_max);

    simdata(i,:) = [psi r delta delta_c psi_w y x_hat(2:3)' ...
        r_hat b_hat bias_true y_r];

    % Wrapped innovation represented by an equivalent held compass input
    x_hat = Ad*x_hat + Bd*[C*x_hat+e; u_obs];
    xs = Ads*xs + Bds*delta_dot + Eds*b_true;
    xw = Adw*xw + Bdw*wave_noise;
    wave_noise = randn/sqrt(h);  % Unit-power white noise held over next step
end

% ------------------------------------------------------------------------------
%% PLOTS
% ------------------------------------------------------------------------------
plotdata = rad2deg(simdata);    

figure(1); clf;
subplot(311)
plot(t,plotdata(:,6),'Color',[0.7 0.7 0.7]); hold on;
plot(t,plotdata(:,1),'b',t,plotdata(:,8),'r',t,rad2deg(psi_ref)*ones(size(t)),'k--');
hold off; grid on;
ylabel('Heading (deg)'); title(['Passive wave filtering: ' sensor_name]);
legend('Compass (LF + WF)','True LF heading','Estimated LF heading','Setpoint');

subplot(312)
plot(t,plotdata(:,5),'b',t,plotdata(:,7),'r'); grid on;
ylabel('WF heading (deg)'); legend('True','Estimated');
xlim([max(0,T_final-100) T_final]); % Show the last 100 s of wave motion

subplot(313)
plot(t,plotdata(:,2),'b',t,plotdata(:,9),'r');
hold off; 
grid on;
xlabel('Time (s)'); ylabel('Yaw rate (deg/s)'); 
switch sensor_option
    case 1
        legend('True LF yaw rate','Estimated yaw rate');
        title('Yaw rate estimated from yaw angle and rudder input');
    case 2
        legend('True LF yaw rate','Filtered and bias-corrected gyro rate');
        title('Measured yaw rate after wave filtering and bias correction');
end

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','axes'),'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',10)

figure(2); clf;
subplot(211)
plot(t,plotdata(:,4),'r--',t,plotdata(:,3),'b'); 
grid on;
ylabel('Rudder (deg)'); title('Rudder response and estimated bias');
legend('Command','Actual');

subplot(212)
plot(t,plotdata(:,11),'b',t,plotdata(:,10),'r'); 
grid on;
xlabel('Time (s)');
switch sensor_option
    case 1
        ylabel('Rudder bias (deg)');
    case 2
        ylabel('Rate-gyro bias (deg/s)');
end 
legend('True bias','Estimated bias');

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','axes'),'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',10)

% Show the large WF gyro signal separately from the LF tracking result.
if sensor_option == 2
    figure(3); clf;
    subplot(211)
    plot(t,plotdata(:,12),'Color',[0.6 0.6 0.6]); grid on;
    ylabel('Yaw rate (deg/s)');
    title('Raw gyro measurement includes wave-frequency motion and bias');

    subplot(212)
    plot(t,plotdata(:,2),'b',t,plotdata(:,9),'r'); grid on;
    xlabel('Time (s)'); ylabel('Yaw rate (deg/s)');
    title('Notch and low-pass filtering reveal the LF yaw rate');
    legend('True LF yaw rate','Filtered and bias-corrected gyro rate');
    set(findall(gcf,'type','line'),'linewidth',1.5)
    set(findall(gcf,'type','axes'),'FontSize',12)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',10)
end

function sensor_option = selectMeasurements
% Large-font selection dialog for classroom demonstrations.
sensor_option = 0;
d = dialog('Name','Passive wave filter - available measurements', ...
    'Tag','PassiveWavefilterSelection','Position',[100 100 800 440], ...
    'Resize','off','Visible','off');
setappdata(d,'sensor_option',0);
uicontrol(d,'Style','text','Position',[25 375 750 40], ...
    'String','Choose available motion measurements', ...
    'FontSize',16,'FontWeight','bold','HorizontalAlignment','left');
uicontrol(d,'Style','text','Position',[25 245 750 110], ...
    'String',sprintf(['YAW ANGLE ONLY\n' ...
    'The observer estimates yaw rate from yaw-angle measurements\n' ...
    'and the known rudder input. No yaw-rate sensor is required.']), ...
    'FontSize',14,'HorizontalAlignment','left');
uicontrol(d,'Style','text','Position',[25 110 750 115], ...
    'String',sprintf(['YAW ANGLE AND YAW RATE\n' ...
    'A gyro measures yaw rate. A notch filter attenuates wave-frequency\n' ...
    'oscillations, and low-pass filtering smooths the signal.\n' ...
    'The observer also estimates and corrects gyro bias.']), ...
    'FontSize',14,'HorizontalAlignment','left');
uicontrol(d,'Style','pushbutton','Position',[25 25 360 65], ...
    'String','Yaw angle only','FontSize',16, ...
    'Callback',@(~,~)choose(1));
uicontrol(d,'Style','pushbutton','Position',[415 25 360 65], ...
    'String','Yaw angle and yaw rate','FontSize',16, ...
    'Callback',@(~,~)choose(2));
movegui(d,'center');
set(d,'Visible','on');
uiwait(d);
if isgraphics(d)
    sensor_option = getappdata(d,'sensor_option');
    delete(d);
end

    function choose(value)
        setappdata(d,'sensor_option',value);
        uiresume(d);
    end
end
