% Passive wave filtering and nonlinear MIMO PID control for dynamic positioning.
% Full-scale 3-DOF vessel using Froude-scaled CyberShip II mass and damping.
% GNSS positions and compass heading contain LF motion + WF motion:
%
%   y = eta + eta_w,    eta = [north east psi]',    nu = [u v r]'.
%
% The 15-state passive observer follows Chapter 13 (Fossen, 2027):
%
%   xi_hat_dot  = Aw*xi_hat + K1*e
%   eta_hat_dot = R*nu_hat + K2*e
%   b_hat_dot   = -inv(T_b)*b_hat + K3*e
%   M*nu_hat_dot = -D*nu_hat + R'*b_hat + tau + R'*K4*e
%
% Here e = y - eta_hat - Cw*xi_hat, with only e(3) wrapped by ssa.
% R is evaluated using the measured heading.
% The environmental bias b is a force/moment vector in NED, not a rudder bias.
% The controller uses estimated LF position, heading and body-fixed velocity.
%
% PIDnonlinearMIMO uses diagonal pole-placement gains and an internal
% reference filter. A 10 deg heading command is applied at 30 s, while
% north and east positions are held at zero. Constant environmental loads
% demonstrate bias estimation during the 600 s simulation.
% Passive observer gains satisfy 1/T_b << K3/K4 < w_0 < w_c in each DOF.
% Only base MATLAB and MSS PIDnonlinearMIMO, sat, ssa, Rzyx and rk4 are used.
%
% All simulation inputs, states, observer gains and time are in full-scale SI
% units. CyberShip II's mass and damping matrices are converted once using
% a 1:70 length ratio: a 1.255 m model represents an 87.85 m ship. With
% G = diag([1 1 scale_ratio]) and rho_ratio = rho_full/rho_model:
%   M = rho_ratio*scale_ratio^3 * G*M_model*G
%   D = rho_ratio*scale_ratio^(5/2) * G*D_model*G
% Thus mass, sway-yaw coupling and yaw inertia scale as s^3, s^4 and s^5;
% linear damping terms scale as s^(5/2), s^(7/2) and s^(9/2).
% simdata contains full-scale SI data; plots convert angles to deg and loads
% to kN/kNm. There is no additional time or motion scaling in the plots.
%
% References:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%      Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%   Ihle, I.-A. F., Skjetne, R. and Fossen, T. I. (2004). Nonlinear Formation
%      Control of Marine Craft with Experimental Results. IEEE Conference on
%      Decision and Control (CDC), Sec. III-A: CyberShip II is a 1:70 scale
%      offshore supply vessel, with model length 1.255 m.
%
% Author:    Thor I. Fossen
% Date:      2026-10-25

clearvars; close all;
clear PIDnonlinearMIMO;         % Reset integral and reference-filter states
rng(23341);

% ------------------------------------------------------------------------------
%% USER INPUTS
% ------------------------------------------------------------------------------
% Full-scale vessel properties obtained from the CyberShip II scale model
L_model = 1.255;                % CyberShip II overall length (m)
scale_ratio = 70;               % Full-scale/model-scale length ratio
L_full = scale_ratio*L_model;   % Equivalent full-scale overall length (87.85 m)
rho_ratio = 1;                  % Full-scale/model-scale fluid density ratio

h = 0.1;                       % Full-scale sampling time (s)
T_final = 1200;                 % Simulation time (s)
t_step = 30;                   % Apply the DP setpoint at this time (s)
eta_ref = [0; 0; deg2rad(10)]; % North (m), east (m), heading (rad)
b_true = [30e3; -15e3; 50e3];  % Constant NED loads [N N Nm]'

% Full-scale vessel: 
%    eta_dot = R(psi)*nu,
%    M*nu_dot + D*nu = tau + R(psi)'*b_true
M_model = [25.8 0 0; 0 33.8 1.0115; 0 1.0115 2.76];
D_model = [2 0 0; 0 7 0.1; 0 0.1 0.5];
G = diag([1 1 scale_ratio]);
M = rho_ratio*scale_ratio^3 * G*M_model*G;
D = rho_ratio*scale_ratio^(5/2) * G*D_model*G;
tau_max = rho_ratio*[scale_ratio^3; scale_ratio^3; scale_ratio^4] ...
    .* [2; 1.5; 1.5];          % Scaled force/moment limits [N N Nm]'

% MIMO PID controller, including its reference filter
wn = 0.05; zeta = 0.9;         % Full-scale closed-loop tuning
T_f = 10;                      % Setpoint low-pass-filter time constant (s)

% Full-scale wave-induced NED position and heading
w_0 = 0.8;                     % Peak frequency (rad/s), period 7.85 s
lambda = 0.1; sigma = 0.5;
wave_gain = [2; 2; deg2rad(3)];
I3 = eye(3); Z3 = zeros(3);
Aw = [Z3 I3; -w_0^2*I3 -2*lambda*w_0*I3];
Bw = [zeros(3,1); 2*lambda*w_0*sigma*wave_gain];
Cw = [Z3 I3];
PHI = expm(h * [Aw Bw; zeros(1,7)]);
Adw = PHI(1:6,1:6); Bdw = PHI(1:6,7);

% Passive observer gains, with finite bias leakage
T_b = 10000*I3; w_c = 1.3*w_0; zeta_n = 1;
K1 = [-2*(zeta_n-lambda)*(w_c/w_0)*I3; ...
       2*w_0*(zeta_n-lambda)*I3];
K2 = w_c*I3;
w_v = 0.05;                   % LF velocity correction bandwidth (rad/s)
K4 = w_c*w_v*diag(diag(M));   % Dimensional force/moment innovation gains
K3 = 0.01*K4;                 % Bias adaptation; K3/K4 = 0.01 1/s
assert(all(1./diag(T_b) < diag(K3)./diag(K4)) && ...
    all(diag(K3)./diag(K4) < w_0) && w_0 < w_c, ...
    'The passive observer gain constraint must be satisfied.');

% Combined plant and observer states:
% x = [eta(3); nu(3); xi_hat(6); eta_hat(3); b_hat(3); nu_hat(3)].
x = zeros(21,1);
xw = zeros(6,1);
t = (0:h:T_final)';
nTimeSteps = length(t);
% Columns: eta, nu, eta_w, y, eta_hat, nu_hat, eta_w_hat, b_hat,
%          applied tau, commanded tau, commanded eta (three columns each).
simdata = zeros(nTimeSteps,33);

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
for i = 1:nTimeSteps

    eta = x(1:3); nu = x(4:6);
    xi_hat = x(7:12); eta_hat = x(13:15);
    b_hat = x(16:18); nu_hat = x(19:21);
    eta_w = Cw*xw;
    y = eta + eta_w;            % GNSS positions and compass heading

    eta_command = eta_ref * (t(i) >= t_step);
    tau_command = PIDnonlinearMIMO(eta_hat,nu_hat,eta_command,M,wn,zeta,T_f,h);
    tau = sat(tau_command, tau_max);

    simdata(i,:) = [eta; nu; eta_w; y; eta_hat; nu_hat; ...
        Cw*xi_hat; b_hat; tau; tau_command; eta_command]';

    % Apply the same saturated force to the vessel and observer
    x = rk4(@vesselDynamics,h,x,tau,b_true,y,M,D,Aw,Cw,T_b,K1,K2,K3,K4);
    xw = Adw*xw + Bdw*randn/sqrt(h);
end

% ------------------------------------------------------------------------------
%% PLOTS
% ------------------------------------------------------------------------------
plotdata = simdata;
angle_columns = [3 6 9 12 15 18 21 33];
plotdata(:,angle_columns) = rad2deg(simdata(:,angle_columns));
plotdata(:,22:30) = plotdata(:,22:30)/1000; % Forces in kN, moments in kNm
position_name = {'North (m)','East (m)','Heading (deg)'};
velocity_name = {'Surge u (m/s)','Sway v (m/s)','Yaw rate (deg/s)'};
load_name = {'Surge force (kN)','Sway force (kN)','Yaw moment (kNm)'};
bias_name = {'North load (kN)','East load (kN)','Yaw moment (kNm)'};

figure(1); clf;
for j = 1:3
    subplot(3,1,j)
    plot(t,plotdata(:,9+j),'Color',[0.7 0.7 0.7]); hold on;
    plot(t,plotdata(:,j),'b',t,plotdata(:,12+j),'r', ...
        t,plotdata(:,30+j),'k--'); hold off; grid on;
    ylabel(position_name{j});
    if j == 1
        title(sprintf('Full-scale DP wave filtering: L = %.2f m',L_full));
        legend('Measured (LF + WF)','True LF','Estimated LF','Setpoint');
    end
end
xlabel('Full-scale time (s)');

figure(2); clf;
for j = 1:3
    subplot(3,1,j)
    plot(t,plotdata(:,6+j),'b',t,plotdata(:,18+j),'r'); grid on;
    ylabel(position_name{j}); xlim([max(0,T_final-100) T_final]);
    if j == 1
        title('Wave-frequency position and heading');
        legend('True WF','Estimated WF');
    end
end
xlabel('Full-scale time (s)');

figure(3); clf;
for j = 1:3
    subplot(3,1,j)
    plot(t,plotdata(:,3+j),'b',t,plotdata(:,15+j),'r'); grid on;
    ylabel(velocity_name{j});
    if j == 1
        title('Body-fixed velocities estimated from position and heading');
        legend('True LF','Estimated LF');
    end
end
xlabel('Full-scale time (s)');

figure(4); clf;
for j = 1:3
    subplot(3,1,j)
    plot(t,plotdata(:,27+j),'r--',t,plotdata(:,24+j),'b'); grid on;
    ylabel(load_name{j});
    if j == 1
        title('DP control forces and moment'); legend('Command','Applied');
    end
end
xlabel('Full-scale time (s)');

figure(5); clf;
for j = 1:3
    subplot(3,1,j)
    plot(t,b_true(j)/1000*ones(size(t)),'b',t,plotdata(:,21+j),'r'); grid on;
    ylabel(bias_name{j});
    if j == 1
        title('Environmental loads estimated in NED');
        legend('True bias','Estimated bias');
    end
end
xlabel('Full-scale time (s)');

figure(6); clf;
plot(plotdata(:,11),plotdata(:,10),'Color',[0.7 0.7 0.7]); hold on;
plot(plotdata(:,2),plotdata(:,1),'b', ...
    plotdata(:,14),plotdata(:,13),'r',eta_ref(2),eta_ref(1),'kx', ...
    'MarkerSize',10); hold off; grid on; axis equal;
xlabel('East (m)'); ylabel('North (m)');
title(sprintf('Full-scale dynamic positioning: L = %.2f m',L_full));
legend('GNSS (LF + WF)','True LF','Estimated LF','Setpoint');

for f = 1:6
    figure(f);
    set(gcf,'Name',sprintf('Full-scale ship: L = %.2f m',L_full), ...
        'NumberTitle','off');
    if any(f == [1 3 4 5])
        set(findall(gcf,'type','axes'),'XLim',[t(1) t(end)]);
    end
    set(findall(gcf,'type','line'),'linewidth',1.5)
    set(findall(gcf,'type','axes'),'FontSize',12)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',10)
end

% ------------------------------------------------------------------------------
% vesselDynamics(x,tau,b_true,y,M,D,Aw,Cw,T_b,K1,K2,K3,K4)
% ------------------------------------------------------------------------------
function xdot = vesselDynamics(x,tau,b_true,y,M,D,Aw,Cw,T_b,K1,K2,K3,K4)
% Full-scale vessel and the book's 15-state nonlinear passive DP observer.
eta = x(1:3); nu = x(4:6);
xi_hat = x(7:12); eta_hat = x(13:15);
b_hat = x(16:18); nu_hat = x(19:21);
R = Rzyx(0,0,eta(3));           % Actual vessel heading (plant only)
R_y = Rzyx(0,0,y(3));           % Measured heading (observer)
e = y - eta_hat - Cw*xi_hat;
e(3) = ssa(e(3));
xdot = [R*nu; M\(-D*nu + tau + R'*b_true); ...
    Aw*xi_hat + K1*e; R_y*nu_hat + K2*e; ...
    -(T_b\b_hat) + K3*e; ...
    M\(-D*nu_hat + R_y'*b_hat + tau + R_y'*K4*e)];
end
