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

clear PIDnonlinearMIMO    % clear the persistent PID variables
clearvars;

t_f = 1000;                     % final simulation time (s)
h   = 0.05;                     % sample time (s)
N = round(t_f/h);               % number of samples

% DP control law
psi_ref = deg2rad(10);
eta_ref = [ 0 0 psi_ref ]';     % DP setpoint
wn = diag([0.02, 0.02, 0.5]);     % closed-loop natural frequencies
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

% Thrust_max(i) = K(i) * n_max(i)^2
n_max = [250, 250, 150, 250, 160, 160]'; % RPM saturation limits

% Tunnel thruster: 3.2 * 250^2 = 200 kN
% Main propeller: 31.2 * 160^2 = 799 kN
K_thr = diag([3.2, 3.2, 3.2, 3.2, 31.2, 31.2]);

T_thr = thrConfig( {'T','T','T','T','M','M'}, ...
    [30, 22, -22, -30, -L/2, -L/2], [0, 0, 0, 0, 5, -5] );

W = diag([1 1 1 1 1 1]);    % Control allocation thruster weights

% Initial states
eta = zeros(3,1);     % eta = [ x y psi ]'
nu = zeros(3,1);      % nu  = [ u v r ]'

%% MAIN LOOP
simdata = zeros(N+1,7);              % memory allocation

for i=1:N+1

    time = (i-1) * h;                % simulation time in seconds

    % Measurements
    eta(1) = eta(1) + 0.01 * randn;
    eta(2) = eta(2) + 0.01 * randn;
    eta(3) = eta(3) + 0.001 * randn;

    % DP control system
    tau_PID = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h);  
    
    % Control allocation using the weight matrix W
    u_alloc = allocPseudoinverse(K_thr,T_thr,W,tau_PID);
    n = sign(u_alloc) .* sqrt(abs(u_alloc));

    % Propeller dynamics and saturation
    for j = 1:length(n)                     % RPM saturation
        if abs(n(j)) > n_max(j)
            n(j) = sign(n(j)) * n_max(j);
        end
    end
    
    u_thr = abs(n) .* n;                     % quadratic RPM

    % Ship dynamics
    tau = T_thr * K_thr * u_thr;
    xdot = supply([eta; nu],tau);  
    R = Rzyx(0,0,eta(3));
    
    % Store data for presentation
    simdata(i,:) = [time,nu',eta']; 
    
   % Euler's integration methods (k+1), (Fossen 2021, Eq. B27-B28)
   % x = x + h * xdot is replaced by forward and backward Euler integration
   nu = nu + h * xdot(4:6);          % Forward Euler 
   eta = eta + h * R * nu;           % Backward Euler

end

%% PLOTS
t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));

figure(1); 
plot(y,x)
title('North-East plot (m)')
xlabel('E'); ylabel('N'); grid
axis equal;
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(221)
plot(t,r)
xlabel('time (s)')
title('yaw rate r (deg/s)')
grid
subplot(222)
plot(t,psi,[0,t(end)],rad2deg([psi_ref psi_ref]))
xlabel('time (s)')
title('yaw angle \psi (deg)')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)


