% Discrete-time Kalman wave filtering and PID heading control of the Mariner 
% class container ship. The compass measures LF heading + WF heading. 
% A five-state Kalman filter separates these components before feedback to 
% the MSS headingAutopilot.
%
% The second-order Nomoto model is
%
%   r(s)/delta(s) = K * (1 + T3*s) / ((1 + T1*s)*(1 + T2*s))
%
% The filter uses T = T1+T2-T3 and x_hat = [xi_w psi_w psi r b]', where
% xi_w is the wave oscillator state and b is an equivalent rudder bias
% modeled as a random walk (b_dot = w_b), with rudder input delta + b:
%
%   x[k+1] = Ad*x[k] + Bd*delta[k] + w[k]
%   y[k]   = C*x[k] + v[k]
%
% This example calls headingAutopilot with its internal reference model and 
% PID tuning. No separate refModel call is needed. Its reference
% bandwidth is set independently using the optional wn_d input. A constant
% 1 deg equivalent rudder bias is added to test disturbance estimation.
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
h = 0.1;                         % Sampling time (s)
T_final = 600;                   % Simulation time (s)
psi_ref = deg2rad(30.0);         % Commanded heading (rad)
b_true = deg2rad(1.0);           % Unknown equivalent rudder bias (rad)

% Mariner ship and rudder servo
K = 0.185;
T1 = 118; T2 = 7.8; T3 = 18.5;  % Second-order Nomoto time constants (s)
T = T1 + T2 - T3;               % First-order Nomoto time constant (s)
T_delta = 1;                    % Rudder servo time constant (s)
delta_max = deg2rad(20);        % Rudder angle limit (rad)
rate_max = deg2rad(5);          % Rudder rate limit (rad/s)

% MSS PID heading autopilot and its internal reference model
wn = 0.1; zeta = 0.5;          % Td = 10 s and Ti = 100 s, as in the case study
wn_d = 0.05;                   % Reference-model bandwidth (rad/s)
r_max = deg2rad(1);            % Desired yaw-rate limit (rad/s)
M_nomoto = T/K;                % Nomoto inertia for a rudder-angle input

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

% Five-state Kalman filter and exact zero-order-hold discretization
A = [0       1               0  0     0
    -w_0^2  -2*lambda*w_0    0  0     0
     0       0               0  1     0
     0       0               0 -1/T   K/T
     0       0               0  0     0];
B = [0 0 0 K/T 0]';
C = [0 1 1 0 0];
PHI = expm(h * [A B; zeros(1,6)]);
Ad = PHI(1:5,1:5);
Bd = PHI(1:5,6);

% Continuous noise intensities. Small yaw-rate process noise makes a steady
% disturbance appear as rudder bias rather than unexplained rate forcing.
% Van Loan's formula gives Qd = integral_0^h exp(A*t)*Qc*exp(A'*t) dt.
% Rd = Rc/h interprets the compass sample as an interval average. These
% are filter tuning covariances; no extra compass noise is added here.
Qc = diag([0 100 0 1e-6 0.01]); Rc = 1;
V = expm(h * [-A Qc; zeros(5) A']);
Qd = Ad * V(1:5,6:10); Qd = (Qd + Qd')/2;
Rd = Rc/h;
x_prd = zeros(5,1); P_prd = eye(5);
I5 = eye(5);

% Initial states and simulation data
xs = zeros(4,1);               % [psi r r_dot delta]'
xw = zeros(2,1);               % [xi_w psi_w]'
t = (0:h:T_final)';
nTimeSteps = length(t);
% Columns: psi, r, delta, delta_c, psi_w, y, psi_w_hat, psi_hat, r_hat,
%          b_hat, b_true. Positive b_true requires negative counter-rudder.
simdata = zeros(nTimeSteps,11);
wave_noise = randn/sqrt(h);

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
for i = 1:nTimeSteps

    % Compass: LF heading + WF heading
    psi = xs(1); r = xs(2); delta = xs(4);
    psi_w = xw(2);
    y = psi + psi_w;

    % Discrete-time Kalman corrector at t(i), with angle wrapping
    Kf = P_prd*C' / (C*P_prd*C' + Rd);
    x_hat = x_prd + Kf * ssa(y - C*x_prd);
    IKC = I5 - Kf*C;
    P_hat = IKC*P_prd*IKC' + Kf*Rd*Kf';  % Joseph covariance update

    % PID heading control uses the estimated LF heading and yaw rate
    delta_c = headingAutopilot(x_hat(3),x_hat(4),psi_ref, ...
        r_max,M_nomoto,wn,zeta,h,wn_d);

    % Rudder servo with angle and rate limits
    delta_ref = sat(delta_c, delta_max);
    delta_dot = sat((delta_ref - delta) / T_delta, rate_max);

    % Store signals at t(i), before advancing to t(i+1)
    simdata(i,:) = [psi r delta delta_c psi_w y x_hat(2:5)' b_true];

    % Kalman predictor: use the measured rudder, not its commanded value
    x_prd = Ad*x_hat + Bd*delta;
    P_prd = Ad*P_hat*Ad' + Qd;
    P_prd = (P_prd + P_prd')/2;

    % Discrete-time ship and wave model propagation
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
ylabel('Heading (deg)'); title('Discrete-time Kalman wave filtering and PID control');
legend('Compass (LF + WF)','True LF heading','Estimated LF heading','Setpoint');

subplot(312)
plot(t,plotdata(:,5),'b',t,plotdata(:,7),'r'); grid on;
ylabel('WF heading (deg)'); legend('True','Estimated');
xlim([max(0,T_final-100) T_final]); % Show the last 100 s of wave motion

subplot(313)
plot(t,plotdata(:,2),'b',t,plotdata(:,9),'r'); 
grid on;
xlabel('Time (s)'); ylabel('Yaw rate (deg/s)'); 
legend('True LF rate','Estimated LF rate');

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
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
xlabel('Time (s)'); ylabel('Rudder bias (deg)'); 
legend('True bias','Estimated bias');

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',10)
