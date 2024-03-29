% SIMsupply  User editable script for simulation of the linear supply
% vessel model, [xdot, U] = supply(x,tau), which returns the speed and the 
% time derivative xdot = A * x + B * tau of the state vector: 
% x = [ x y psi u v r ]' for a supply vessel length L = 76 m.
%
% Calls:      supply.m
%             PIDnonlinearMIMO.m
%
% Author:     Thor I. Fossen
% Date:       2024-03-28
% Revisions:  

clear PIDnonlinearMIMO          % clear the persistent PID variables
clearvars;

t_f = 500;                      % final simulation time (s)
h   = 0.05;                     % sample time (s)
N = round(t_f/h);               % number of samples

% Environment
Vc = 0.5;                       % current speed (m/s)
betaVc = deg2rad(20);           % current direction (rad)

% DP control law
psi_ref = deg2rad(50);
eta_ref = [ 5 5 psi_ref ]';     % DP setpoint
wn = diag([0.1, 0.1, 0.2]);     % closed-loop natural frequencies
zeta = diag([1, 1, 1]);         % closed-loop relative damping ratios
T_f = 10;                       % setpoint LP-filter time constant

M = 1e9 * [...                  % supply vessel mass matrix
    0.0068         0         0
         0    0.0113   -0.0340
         0   -0.0340    4.4524 ];

% Thrust coefficient and configuration matrices (Fossen 2021, Ch. 11.2)
%   #1 Bow tunnel thruster (RPM)
%   #2 Bow tunnel thruster (RPM)
%   #3 Stern tunnel thruster (RPM)
%   #4 Stern tunnel thruster (RPM)
%   #5 Right main propeller (RPM)
%   #6 Left main propeller (RPM)
L    =  76.2;                           % length of ship (m)
T_n = 1;                                % propeler speed time constant (s)

% Thrust_max(i) = K(i) * n_max(i)^2
n_max = [250, 250, 250, 250, 160, 160]'; % RPM saturation limits

% Tunnel thruster: 3.2 * 250^2 = 200 kN
% Main propeller: 31.2 * 160^2 = 799 kN
K_thr = diag([3.2, 3.2, 3.2, 3.2, 31.2, 31.2]);

T_thr = thrConfig( {'T','T','T','T','M','M'}, ...
    [30, 22, -22, -30, -L/2, -L/2], [0, 0, 0, 0, 8, -8] );

W = diag([1 1 1 1 1 1]);    % Control allocation propeller weights

% Initial states
eta = zeros(3,1);     % eta = [ x y psi ]'
nu = zeros(3,1);      % nu  = [ u v r ]'
n = zeros(6,1);       % vector of propeller speed states

%% MAIN LOOP
simdata = zeros(N+1,19);              % memory allocation

for i=1:N+1

    time = (i-1) * h;                % simulation time in seconds

    % Measurements
    eta(1) = eta(1) + 0.01 * randn;
    eta(2) = eta(2) + 0.01 * randn;
    eta(3) = eta(3) + 0.0001 * randn;

    % Current velocities
    u_c = Vc * cos(betaVc - eta(3));    % current surge velocity
    v_c = Vc * sin(betaVc - eta(3));    % current sway velocity
    nu_r = nu - [u_c, v_c, 0]';         % relative velocity vector

    % DP control system
    tau_PID = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h);

    % Control allocation using the weight matrix W
    u_alloc = allocPseudoinverse(K_thr,T_thr,W,tau_PID);
    n_c = sign(u_alloc) .* sqrt(abs(u_alloc));

    % Propeller saturation
    for j = 1:6                     % RPM saturation
        if abs(n(j)) > n_max(j)
            n(j) = sign(n(j)) * n_max(j);
        end
    end

    % Ship and propeller dynamics 
    u_thr = abs(n) .* n;
    tau = T_thr * K_thr * u_thr;
    xdot = supply([eta; nu_r], tau);
    ndot = (n_c - n) / T_n;

    % Store data for presentation
    simdata(i,:) = [time,nu',eta',n',n_c'];

    % Euler's integration methods (k+1), (Fossen 2021, Eq. B27-B28)
    % x = x + h * xdot is replaced by forward and backward Euler integration
    nu = nu + h * xdot(4:6);                % Forward Euler
    eta = eta + h * Rzyx(0,0,eta(3)) * nu;  % Backward Euler
    n = n + h * ndot;

end

%% PLOTS
t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));
n     = simdata(:,8:13);
n_c   = simdata(:,14:19);

figure(1); 
plot(y,x)
title('North-East plot (m)')
xlabel('E'); ylabel('N'); grid
axis equal;
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(211)
plot(t,r)
xlabel('time (s)')
title('yaw rate r (deg/s)'); grid
subplot(212)
plot(t,psi,[0,t(end)],rad2deg([psi_ref psi_ref]))
xlabel('time (s)')
title('yaw angle \psi (deg)'); grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(3)
subplot(321)
plot(t,n(:,1),t,n_c(:,1),[0 t(end)],...
    [n_max(1) n_max(1)],'k',[0 t(end)],-[n_max(1) n_max(1)],'k')
xlabel('time (s)')
title('Bow tunnel thruster n_1 (RPM)'); grid
subplot(322)
plot(t,n(:,2),t,n_c(:,2),[0 t(end)],...
    [n_max(2) n_max(2)],'k',[0 t(end)],-[n_max(2) n_max(2)],'k')
xlabel('time (s)')
title('Bow tunnel thruster n_2 (RPM)'); grid
subplot(323)
plot(t,n(:,3),t,n_c(:,3),[0 t(end)],...
    [n_max(3) n_max(3)],'k',[0 t(end)],-[n_max(3) n_max(3)],'k')
xlabel('time (s)')
title('Stern tunnel thruster n_3 (RPM)'); grid
subplot(324)
plot(t,n(:,4),t,n_c(:,4),[0 t(end)],...
    [n_max(4) n_max(4)],'k',[0 t(end)],-[n_max(4) n_max(4)],'k')
xlabel('time (s)')
title('Stern tunnel thruster n_4 (RPM)'); grid
subplot(325)
plot(t,n(:,5),t,n_c(:,5),[0 t(end)],...
    [n_max(5) n_max(5)],'k',[0 t(end)],-[n_max(5) n_max(5)],'k')
xlabel('time (s)')
title('Main propeller n_5 (RPM)'); grid
subplot(326)
plot(t,n(:,6),t,n_c(:,6),[0 t(end)],...
    [n_max(6) n_max(6)],'k',[0 t(end)],-[n_max(6) n_max(6)],'k')
xlabel('time (s)')
title('Main propeller n_6 (RPM)'); grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)


