% SIMosv 
% User editable script for simulation of an Offshore Supply Vessel ('osv.m'). 
% The OSV is expressed by the nonlinear equations of motion based on
% Fossen (2021, Eqs. 6.111-6.116):
%   
%   eta_dot = J(eta) * nu
%   nu_dot = nu_c_dot + Minv * (tau_thr +  tau_drag + tau_crossflow...
%          - (CRB + CA + D) * nu_r - G * eta )
%
% where Minv = inv(MRB + MA) and nu_r = nu - nu_c is the relative velocity
% vector. The generalized thrust, tau_thr = T(alpha) * K_thr * u_thr, is 
% defined by 
%
%   alpha = [alpha1, alpha2]'  - azimuth angles
%   u_thr = [bowThruster 1, bowThruster 2, sternAzimuth 1, sternAzimuth 2]
%
% where u_thr(i) = abs(n(i)) * n(i) for i = 1,...,4 is the squared propeller 
% speed. The DP control law is chosen as a MIMO nonlinear PID controller 
% for setpoint regulation using the function PIDnonlinearMIMO.m based on 
% Fossen (2021, Algorithm 15.2) where
%   
%   z_int = z_int + h * (eta - eta_d) 
%   tau_thr = -R(psi)' * ( Kp * (eta - eta_d) + Kd * nu + Ki * z_int )
% 
% Both unconstrained control alloction (pseduoinverse) and constrained
% control allocation (dynamic optimization) are implemented using Fossen
% (2021, Sections 11.2.2-11.2.3).
%
% Calls: 
%   fmincon.m            - Sequential Quadratic Programming (SQP).
%                          Requires the MATLAB optimization toolbox.
%   osv.m                - OSV equations of motion.
%   PIDnonlinearMIMO.m   - MIMO nonlinear PID controller.
%   allocPseudoinverse.m - Unconstrained control allocation. 
%   optimalAlloc.m       - Custom function for constrained control 
%                          allocation, defined within the 'SIMosv.m' script.
%
% It is possible to run SIMosv without the MATLAB optimization toolbox 
% by choosing ALLOC = 0. Then the control allocation problem is solved with
% two constant azimuth angles: alpha0 = deg2rad([-28 28])', which
% corresponds to the minimum condition number of the thruster configuration
% matrix T(alpha). The constrained control allocation problem (ALLOC = 1) 
% is a nonlinear dynamic optimization problem, which requires the Matlab 
% optimization toolbox.
%
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:

clear PIDnonlinearMIMO;      % clear persitent variables
clearvars;                   % clear other variables
osv;                         % display the OSV main data

% Control allocation 
ALLOC = 1;                   % 0 = constant azimuth, 1 = dynamic optimization
alpha0 = deg2rad([-28 28])'; % constant azimuth, minimizing cond(T_thr)

% Simulation parameters
h = 0.1;                     % sampling time (s)
N = 2000;                    % number of samples

% DP setpoints at initial time     
x_ref = 0;                                 % North position (m)
y_ref = 0;                                 % East position (m)
psi_ref = deg2rad(0);                      % yaw angle (rad)
eta_ref = [x_ref y_ref psi_ref]'; 

% Ship model parameters
L = 83;                      % length (m)
B = 18;                      % beam (m)
T = 5;                       % draft (m)

% Ocean current
Vc = 1;                      % current speed (m/s)
betaVc = deg2rad(-140);      % current direction (rad)

% Thrust: T = K_max * abs(n/n_max) * (n/n_max) = K_thr * abs(n) * n
K_max = diag([300e3 300e3 655e3 655e3]);   % max propeller thrust (N)
n_max = [140 140 150 150]';                % max propeller speed (rpm)
K_thr = K_max./n_max.^2;                   % thruster coefficient matrix
l_x = [37, 35, -L/2, -L/2];                % thruster x-coordinates
l_y = [0, 0, 7, -7];                       % thruster y-coordinates

% Unconstrained control allocation: Constant azimuth angles, alpha_0   
T_thr = thrConfig( {'T', 'T', alpha0(1),alpha0(2)}, l_x, l_y);

% Constrained control allocation: Dynamic optimization
az_max = deg2rad(60);  % maximum azimuth angles

% Lower and upper bounds: x = [az1 az2 u1 u2 u3 u4 s1 s2 s3] 
% 2 azimuths, 4 quadratic propeller speeds: u = abs(n/n_max) * (n/n_max), 
% and 3 slack variables
lb = [-az_max -az_max -1 -1 -1 -1 -inf -inf -inf];
ub = [ az_max  az_max  1  1  1  1  inf  inf  inf];

alpha_old = alpha0;  % initial vales for dynamic optimization
u_old = [0 0 0 0]';

% MIMO nonlinear PID controller:
% tau = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h)
M = 1e9 * ...  % computed in osv.m
    [ 0.0060         0         0         0   -0.0060         0
          0    0.0080         0    0.0100         0   -0.0284
          0         0    0.0130         0    0.2554         0
          0    0.0100         0    0.3067         0   -0.1299
    -0.0060         0    0.2554         0    6.4508         0
          0   -0.0284         0   -0.1299         0    3.3996 ];
wn = 0.1 * diag([1 1 3]);         % closed-loop natural frequencies (rad/s)
zeta = 1.0 * diag([1 1 1]);       % closed-loop relative damping ratios (-)
T_f = 50;                         % LP-filter time constant (s)

% Initial states
eta = [0 0 0 deg2rad(5) deg2rad(1.3) 0]';  % eta = [x y z phi theta psi]' 
nu  = [0 0 0 0 0 0]';                      % nu  = [u v w p q r]'

% Allocate empty table for simulation data
simdata = zeros(N+1,19); 

% Create a waitbar
h_waitbar = waitbar(0, 'Processing...');
tic;

%% MAIN LOOP
for i = 1:N+1
   
   t = (i-1) * h;                   % time (s)  

   % Update the waitbar
    if mod(i, 10) == 0
        waitbar(i/N, h_waitbar, sprintf('Progress: %3.0f%%', i/N*100));
    end

   % Measurements with white-noise
   eta(1) = eta(1) + 0.06 * randn;     
   eta(2) = eta(2) + 0.06 * randn;
   eta(6) = eta(6) + 0.001 * randn;    

   % MIMO nonlinear PID controller
   if t > 50 
       eta_ref = [0 0 deg2rad(40)]'; 
   end

   tau = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h); 

   % Thrust: T = K_max * abs(n/n_max) * (n/n_max) = K_thr * abs(n) * n
   % Normalized squared propeller speed: u_c = abs(n_c/n_max) * (n_c/n_max)
   if ALLOC == 0    % Unconstrained control allocation: pseduoinverse
   
       alpha_c = alpha0;  % constant azimuth angles
       u_c = allocPseudoinverse( K_max,T_thr,eye(4),tau([1:2,6]) );

   else  % Constrained control allocation: dynamic optimization 

       [alpha_c, u_c, slack] = optimalAlloc(tau([1:2,6]), lb, ub, ...
           alpha_old, u_old, l_x, l_y, K_max, h);
       alpha_old = alpha_c;
       u_old = u_c;

   end

   % Controls: ui = [ n_c(1) n_c(2) n_c(3) n_c(4) alpha_c(1) alpha_c(2) ]'
   u_c = n_max.^2 .* u_c;                   
   n_c = sign(u_c) .* sqrt( abs(u_c) );
   ui = [n_c; alpha_c];

   % Store simulation data in a table   
   simdata(i,:) = [t eta' nu' ui']; 

   % Offshore supply vessel dynamics
   xdot = osv([nu; eta],ui,Vc,betaVc);
   Jmtrx = eulerang(eta(4),eta(5),eta(6));

   % Propagate the vehicle dynamics (k+1), (Fossen 2021, Eq. B27-B28)
   % x = x + h * xdot is replaced by forward and backward Euler integration
   nu = nu + h * xdot(1:6);        % Forward Euler, velocity
   eta = eta + h * Jmtrx * nu;     % Backward Euler, position
   
end

close(h_waitbar);                   % close the waitbar after the loop

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
legend('Yaw angle (deg)','Location','best')
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
legend('p (deg/s)','q (deg/s)','r (deg/s)','Location','best')
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
legend('n_1','n_2','Location','best')

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
legend('n_3','n_4','Location','best')

subplot(313)
hold on;
plot(t,rad2deg(a1),'b','linewidth',2)
plot(t,rad2deg(a2),'r','linewidth',2)
plot([0,t(end)],rad2deg([az_max,az_max]),'k','linewidth',1)
plot([0,t(end)],rad2deg([-az_max,-az_max]),'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Azimuuth angles (deg)'),grid
legend('\alpha_1','\alpha_2','Location','best')

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)


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

% mMximum azimuth angle and pitch rates
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




