function SIMosv()
% SIMosv is compatibel with MATLAB and incorporates dynamic and static 
% optimization techniques for control allocation, though dynamic 
% optimization is not supported in GNU Octave. The script simulates an 
% Offshore Supply Vessel (OSV)  utilizing a Dynamic Positioning (DP) system 
% for stationkeeping and low-speed maneuvering under the influence of ocean 
% currents. The OSV's behavior is modeled by nonlinear equations of motion 
% as specified in Fossen (2021), which includes the following equations:
%
%   eta_dot = J(eta) * nu
%   nu_dot = nu_c_dot + Minv * (tau_thr + tau_drag + tau_crossflow ...
%            - (CRB + CA + D) * nu_r - G * eta)
%
% where:
%   - Minv = inv(MRB + MA) is the inverse of the system mass matrix.
%   - nu_r = nu - nu_c represents the relative velocity vector.
%   - tau_thr = T(alpha) * K_thr * u_thr describes the generalized thrust,
%     with alpha representing azimuth angles and u_thr the propeller speeds.
%
% The DP control strategy employs a MIMO nonlinear PID controller for 
% setpoint regulation, based on Fossen (2021, Algorithm 15.2). The control 
% laws include:
%
%   z_int = z_int + h * (eta - eta_d)
%   tau_thr = -R(psi)' * (Kp * (eta - eta_d) + Kd * nu + Ki * z_int)
%
% Control allocation is implemented both as unconstrained (using 
% pseudoinverse methods) and constrained (via dynamic optimization) 
% techniques, detailed in Fossen (2021, Sections 11.2.2-11.2.3).
%
% Dependencies:
%   This script requires the MATLAB optimization toolbox for dynamic 
%     optimization features due to the use of fmincon for sequential 
%     quadratic programming (SQP). 
%   In Octave, set ALLOC = 0 for static optimization using constant 
%     azimuth angles to minimize the condition number of the thruster 
%     configuration matrix T(alpha).
%
% Main Functions Called:
%   fmincon.m            - For dynamic optimization using SQP.
%   osv.m                - Calculates the OSV's equations of motion.
%   PIDnonlinearMIMO.m   - Implements the MIMO nonlinear PID controller.
%   allocPseudoinverse.m - Performs unconstrained control allocation.
%   optimalAlloc.m       - Conducts constrained control allocation within 
%                          this script.
% Example Usage:
%   Set ALLOC = 0 to run without MATLAB's optimization toolbox.
%   For dynamic optimization, ensure ALLOC = 1 and MATLAB is used.
%
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:
%   None

clear PIDnonlinearMIMO; % Clear persistent variables from previous sessions.
clearvars;              % Clear all other variables.
close all;              % Close all windows.
osv;                    % Initialize or display the OSV's main data.

% Define control allocation mode:
ALLOC = 1;                   % 0 = Use constant azimuths 
                             % 1 = Use dynamic optimization (MATLAB only)

% Constant azimuth angles minimizing the condition number of T_thr                             
alpha0 = deg2rad([-28; 28]); 

% Define simulation parameters:
h = 0.1;                     % Sampling time in seconds.
N = 2000;                    % Number of simulation samples.

% Define DP setpoints:
x_ref = 0;                   % Reference North position in meters.
y_ref = 0;                   % Reference East position in meters.
psi_ref = deg2rad(0);        % Reference yaw angle in radians.
eta_ref = [x_ref, y_ref, psi_ref]';  % Reference positions and heading.

% Vessel and environmental parameters:
L = 83;                      % Vessel length in meters.
B = 18;                      % Vessel beam in meters.
T = 5;                       % Vessel draft in meters.
Vc = 1;                      % Ocean current speed in meters/second.
betaVc = deg2rad(-140);      % Ocean current direction in radians.

% Thruster configuration parameters:
K_max = diag([300e3, 300e3, 655e3, 655e3]);  % Max thrust for each propeller in Newtons.
n_max = [140, 140, 150, 150]';               % Max propeller speeds in RPM.
K_thr = K_max ./ n_max.^2;                   % Thrust coefficient matrix.
l_x = [37, 35, -L/2, -L/2];                  % X-coordinates of thrusters.
l_y = [0, 0, 7, -7];                         % Y-coordinates of thrusters.

% Thruster configuration matrix
T_thr = thrConfig({'T', 'T', alpha0(1), alpha0(2)}, l_x, l_y);  

% Dynamic optimization setup (if applicable):
az_max = deg2rad(60);  % Max azimuth rotation angle in radians.

% Bounds for control variables:
lb = [-az_max, -az_max, -1, -1, -1, -1, -inf, -inf, -inf];  % Lower bounds.
ub = [az_max, az_max, 1, 1, 1, 1, inf, inf, inf];           % Upper bounds.

alpha_old = alpha0;    % Initial values for dynamic optimization.
u_old = [0, 0, 0, 0]'; % Initial propeller speeds.

% Initialize the nonlinear MIMO PID controller:
[~,~,M] = osv();             % OSV 6x6 mass matrix
wn = 0.1 * diag([1 1 3]);    % Natural frequencies for PID tuning.
zeta = 1.0 * diag([1 1 1]);  % Damping ratios for PID tuning.
T_f = 50;                    % Time constant for the low-pass filter in seconds.

% Initialize state vectors for the simulation:
eta = [0, 0, 0, deg2rad(5), deg2rad(1.3), 0]';  % Euler angles and positions.
nu = [0, 0, 0, 0, 0, 0]';                       % Velocity vector.

% Allocate memory for simulation data storage:
simdata = zeros(N+1, 19);  % Pre-allocate matrix for efficiency.

% Octave can only use ALLOC = 0
if isoctave
    ALLOC = 0;
    disp('Control allocation in Octave works only for constant azimuth angles')
end

% Create a progress indicator
h_waitbar = waitbar(0, 'Processing...');    % Display a wait bar 
tic;  % Start a timer to measure the simulation's execution time.

%% Main simulation loop:
for i = 1:N+1
   
   t = (i-1) * h;  % Current simulation time.

   % Update the progress bar every 10 iterations:
   if mod(i, 10) == 0
       waitbar(i/N, h_waitbar, sprintf('Progress: %3.0f%%', i/N*100));
   end

   % Simulate sensor noise and disturbances:
   eta(1) = eta(1) + 0.06 * randn;   % Simulate noise in the North position.
   eta(2) = eta(2) + 0.06 * randn;   % Simulate noise in the East position.
   eta(6) = eta(6) + 0.001 * randn;  % Simulate noise in the yaw angle.

   % Control logic based on the elapsed simulation time:
   if t > 50 
       eta_ref = [0, 0, deg2rad(40)]';  % Change setpoint after 50 seconds.
   end

   % Calculate control forces using the nonlinear MIMO PID controller:
   tau = PIDnonlinearMIMO(eta, nu, eta_ref, M, wn, zeta, T_f, h); 

   % Determine thrust settings based on the selected allocation method:
   if ALLOC == 0  % Unconstrained control allocation using pseudoinverse method.
       alpha_c = alpha0 ; % Use constant azimuth angles.
       u_c = allocPseudoinverse(K_max, T_thr, eye(4), tau([1:2, 6]));  
   else  % Constrained control allocation using dynamic optimization.
       [alpha_c, u_c, ~] = optimalAlloc(tau([1:2, 6]), lb, ub, alpha_old, ...
           u_old, l_x, l_y, K_max, h);
       alpha_old = alpha_c;  % Update for next iteration.
       u_old = u_c;          % Update for next iteration.
   end

   % Controls: ui = [ n_c(1) n_c(2) n_c(3) n_c(4) alpha_c(1) alpha_c(2) ]'
   u_c = n_max.^2 .* u_c;  % Scale control efforts to actual propeller speeds.
   n_c = sign(u_c) .* sqrt(abs(u_c));  % Calculate each propeller's speed.
   ui = [n_c; alpha_c];  

   % Store simulation results:
   simdata(i, :) = [t, eta', nu', ui'];  % Log data for later analysis.

   % Calculate the OSV's dynamics for the next timestep:
   xdot = osv([nu; eta], ui, Vc, betaVc);    
   Jmtrx = eulerang(eta(4), eta(5), eta(6));  

   % Integrate to find the next state:
   nu = nu + h * xdot(1:6);     % Update velocity using forward Euler method.
   eta = eta + h * Jmtrx * nu;  % Update position using backward Euler method.
   
end

close(h_waitbar);  % Close the progress indicator.

%% PLOTS
t     = simdata(:,1);  

x     = simdata(:,2); 
y     = simdata(:,3); 
z     = simdata(:,4); 
phi   = ssa(simdata(:,5)); 
theta = ssa(simdata(:,6)); 
psi   = ssa(simdata(:,7)); 

u     = simdata(:,8);        
v     = simdata(:,9); 
w     = simdata(:,10); 
p     = simdata(:,11);        
q     = simdata(:,12); 
r     = simdata(:,13); 

U = sqrt( u.^2 + v.^2 );        % vessel speed (m/s)

n1 = simdata(:,14);             % propeller speeds (rpm)
n2 = simdata(:,15);            
n3 = simdata(:,16);
n4 = simdata(:,17);

a1 = simdata(:,18);             % propeller azimuth angles (rad)
a2 = simdata(:,19); 

legendLocation = 'best';
if isoctave; legendLocation = 'northeast'; end

%% Position and Euler angle plots
figure(2); clf;
figure(gcf)
subplot(321),plot(y,x)
xlabel('East (m)')
ylabel('North (m)')
title('North-East positions (m)'),grid
subplot(322),plot(t,z)
xlabel('time (s)'),title('Down position (m)'),grid
subplot(312),plot(t,rad2deg(phi),t,rad2deg(theta))
xlabel('time (s)'),title('Roll and pitch angles (deg)'),grid
legend('Roll angle (deg)','Pitch angle (deg)')
subplot(313),plot(t,rad2deg(psi))
xlabel('time (s)'),title('Heading angle (deg)'),grid
legend('Yaw angle (deg)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Velocity plots
figure(3); clf;
figure(gcf)
subplot(311),plot(t,U)
xlabel('time (s)'),title('Speed (m/s)'),grid
subplot(312),plot(t,u,t,v,t,w)
xlabel('time (s)'),title('Linear velocities (m/s)'),grid
legend('u (m/s)','v (m/s)','w (m/s)')
subplot(313),plot(t,rad2deg(p),t,rad2deg(q),t,rad2deg(r))
xlabel('time (s)'),title('Angular velocities (deg/s)'),grid
legend('p (deg/s)','q (deg/s)','r (deg/s)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Propeller plots
figure(4); clf;
figure(gcf)
subplot(311)
hold on;
plot(t,n1,'b','linewidth',2)
plot(t,n2,'r','linewidth',2)
plot([0,t(end)],[n_max(1),n_max(1)],'k','linewidth',1)
plot([0,t(end)],[-n_max(1),-n_max(1)],'k','linewidth',1)
plot([0,t(end)],[n_max(2),n_max(2)],'k','linewidth',1)
plot([0,t(end)],[-n_max(2),-n_max(2)],'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Bow thrusters (rpm)'),grid
legend('n_1','n_2','Location',legendLocation)

subplot(312)
hold on;
plot(t,n3,'b','linewidth',2)
plot(t,n4,'r','linewidth',2)
plot([0,t(end)],[n_max(3),n_max(3)],'k','linewidth',1)
plot([0,t(end)],[-n_max(3),-n_max(3)],'k','linewidth',1)
plot([0,t(end)],[n_max(4),n_max(4)],'k','linewidth',1)
plot([0,t(end)],[-n_max(4),-n_max(4)],'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Stern azimuth thrusters (rpm)'),grid
legend('n_3','n_4','Location',legendLocation)

subplot(313)
hold on;
plot(t,rad2deg(a1),'b','linewidth',2)
plot(t,rad2deg(a2),'r','linewidth',2)
plot([0,t(end)],rad2deg([az_max,az_max]),'k','linewidth',1)
plot([0,t(end)],rad2deg([-az_max,-az_max]),'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Azimuuth angles (deg)'),grid
legend('\alpha_1','\alpha_2','Location',legendLocation)

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

end

%% FUNCTIONS FOR DYNAMIC OPTIMIZATION
function [alpha_opt, u_opt, slack] = ...
    optimalAlloc(tau, lb, ub, alpha_old, u_old, l_x, l_y, K_thr, h)
% [alpha_opt, u_opt, slack] = ...
%    optimalAlloc(tau, lb, ub, alpha_old, u_old, l_x, l_y, K_thr, h)
% Constrained allocation is implemented by using sequential quadratic
% programming (SQP) to find the optimal (alpha, u) values when there
% are amplitude and rate saturations. The mathematical model is
%
% tau = T_thr(alpha) * K_thr * u
%   u = abs(n/n_max) * (n/n_max) - normalized squared propeller speed
%
% where the six controls are
%
% alpha = [alpha1, alpha2]'                        azimuth angles
% u = [bowAzimuth, starboardAzimuth, portAzimuth]  squared propeller speed
%
% Inputs:
%  tau: vector of generalized forces (surge, sway, yaw)
%  lb: vector of lower bounds
%  ub: vector of upper bounds
%  alpha_old: previous values of alpha
%  u_old: previous values of u
%  l_x: vector of x-coordinates, thruster lever arms
%  l_y: vector of y-coordinates, thruster lever arms
%  K_thr: diagonal matrix of maximum thrust
%  h: sampling time
%
% Outputs:
%  alpha_opt: optimal alpha values
%  u_opt: optimal u values, normalized
%  slack: norm of slack variables s = tau - T_thr(alpha) * K_thr * u
%
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:

% Objective function
fun = @(x) objectiveFunction(x,alpha_old,u_old);

% Initial guess: x0 = [az1 az2 u1 u2 u3 u4 s1 s2 s3]
x0 = [deg2rad(-28) deg2rad(28)  0 0 0 0  0 0 0];

% Use fmincon and sequential quadratic programming (SQP) to optimize
nonlcon = @(x) constraints(x, alpha_old, u_old, l_x, l_y, K_thr, tau, h);
options = optimoptions('fmincon','Display','off','Algorithm','sqp');
x_opt = fmincon(fun, x0, [], [], [], [], lb, ub, nonlcon, options);

% Extract results
alpha_opt = x_opt(1:2)';
u_opt = x_opt(3:6)';
slack = norm( x_opt(7:9) );

end

% Cost function
function cost = objectiveFunction(x,alpha_old,u_old)

alpha = x(1:2)';    % azimuth angles
u = x(3:6)';        % quadratic controls
s = x(7:9)';        % slack variables: tau = T(alpha) * K_thr * u + s

w1 = 1;             % weight for squared u
w2 = 1000;          % weight for squared slack variable s
w3 = 10;            % weight for squared change in alpha
w4 = 10;            % weight for squared change in u

cost = w1 * norm(u)^2 ...
    + w2 * norm(s)^2 ...
    + w3 * norm(alpha - alpha_old)^2 ...
    + w4 * norm(u - u_old)^2;

end

% Constraints
function [c, ceq] = constraints(x, alpha_old, u_old, l_x, l_y, K_thr,tau,h)

alpha = x(1:2)';     % azimuth angles
u = x(3:6)';         % quadratic controls
s = x(7:9)';         % slack variables: tau = T(alpha) * K_thr * u + s

% Maximum azimuth angle and pitch rates
max_rate_alpha = 0.3;           % 0.3 rad/s = 17.2 deg/s
max_rate_u = 0.3;               

c(1) = abs( alpha(1) - alpha_old(1) ) / h - max_rate_alpha;
c(2) = abs( alpha(2) - alpha_old(2) ) / h - max_rate_alpha;

c(3) = abs( u(1) - u_old(1) ) / h - max_rate_u;
c(4) = abs( u(2) - u_old(2) ) / h - max_rate_u;
c(5) = abs( u(3) - u_old(3) ) / h - max_rate_u;
c(6) = abs( u(4) - u_old(4) ) / h - max_rate_u;

% Equality constraint
T_alpha = thrConfig( {'T', 'T', alpha(1), alpha(2)}, l_x, l_y);

ceq = T_alpha * K_thr * u - tau + s;

end
