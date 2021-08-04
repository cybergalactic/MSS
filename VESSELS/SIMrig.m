% SIMrig   User editable script for simulation of a semisubmersible with
% two pontoons and four columns. The data is generated using the script
% rig.m, which calls data_rig.m. 
%
%   m =  27 162 500 kg (mass)
%   H_p = 13.0 m (pontoon height)
%   B_p = 12.0 m (pontoon width)
%   L_p = 84.6 m (pontoon length)
%   r_CG = [-0.5 0 -19.5]', CG with respect to the CO
%
% Author:      Thor I. Fossen
% Date:        2019-03-12
% Revisions:   2020-10-22 retuning of PID controller
%              2020-08.04 major update

clearvars

f = 20;              % sampling frequency
h = 1/f;             % sampling time
N = round(1000*f);   % number of samples

% environmnetal data
A_wave = 5000000 * [0.5 0.1 0 100 10 10]';      % wave gain 
nu_c = [-0.3 -0.2 0 0 0 0]';                    % ocean current velocities

% initial values
eta = zeros(6,1);   % generalized position
nu  = zeros(6,1);   % generalized velocity
z3 = zeros(3,1);    % integral state

% Semisub model matrices
MRB = 1.0e+10 * [
    0.0027         0         0         0   -0.0530         0
    0    0.0027         0    0.0530         0   -0.0014
    0         0    0.0027         0    0.0014         0
    0    0.0530         0    3.4775         0   -0.0265
    -0.0530         0    0.0014         0    3.8150         0
    0   -0.0014         0   -0.0265         0    3.7192];

MA = 1.0e+10 * [
    0.0017         0         0         0   -0.0255         0
    0    0.0042         0    0.0365         0         0
    0         0    0.0021         0         0         0
    0    0.0365         0    1.3416         0         0
    -0.0255         0         0         0    2.2267         0
    0         0         0         0         0    3.2049];

D = 1.0e+09 * [
    0.0004         0         0         0   -0.0085         0
    0    0.0003         0    0.0067         0   -0.0002
    0         0    0.0034         0    0.0017         0
    0    0.0067         0    4.8841         0   -0.0034
    -0.0085         0    0.0017         0    7.1383         0
    0   -0.0002         0   -0.0034         0    0.8656];

G = 1.0e+10 * [
    0         0         0         0         0         0
    0         0         0         0         0         0
    0         0    0.0006         0         0         0
    0         0         0    1.4296         0         0
    0         0         0         0    2.6212         0
    0         0         0         0         0         0 ];

M = MRB + MA;    
Minv = inv(M);  

% Nonlinear PID controller  
M3 = diag( [M(1,1), M(2,2), M(6,6) ]);  % 3-DOF matrices 
D3 = diag( [D(1,1), D(2,2), D(6,6) ]);  
G3 = diag( [G(1,1), G(2,2), G(6,6) ]); 

w_n = 0.1;                          % Bandwidth
Kp = M3 * w_n^2  * eye(3) - G3;     % Proportional gain 
Kd = 2 * w_n * M3 - D3;             % Derivative  gain 
Ki = (w_n / 10) * Kp;               % Integral  gain     

%% *************** MAIN SIMULATION LOOP ************************
for i=1:N+1
    t = (i-1)*h;                   % simulation time in seconds

    % kinematics
    phi = eta(4); theta = eta(5); psi = eta(6);
    J = eulerang(phi,theta,psi); 
    R = Rzyx(phi,theta,psi); 
    
    % 3-DOF nonlinear DP control system (no wave filtering)
    eta3 = [eta(1) eta(2) ssa( eta(6) )]' ;
    nu3  = [nu(1) nu(2) nu(6)]';   
    tau3 = -R' * (Kp * eta3 + Ki * z3) - Kd * nu3; 
    tau = [tau3(1) tau3(2) 0 0 0 tau3(3)]';
    
    % equations of motion
    t_wave = A_wave * sin(0.2 * t);
    etadot = J * nu;
    nudot = Minv * (tau - D * (nu-nu_c) - G * eta + t_wave);
    
    % store data for presentation
    xout(i,:) = [t, eta', nu'];
    
    % numerical integration
    eta = euler2(etadot, eta, h);             % Euler's method
    nu  = euler2(nudot, nu, h);     
    z3  = euler2(eta3, z3, h);

end
% *************** END SIMULATION LOOP ************************

% Time-series
t     = xout(:,1);
x     = xout(:,2); 
y     = xout(:,3);          
z     = xout(:,4);   
phi   = xout(:,5) * 180/pi; 
theta = xout(:,6) * 180/pi; 
psi   = xout(:,7) * 180/pi;
u     = xout(:,8);
v     = xout(:,9);
w     = xout(:,10);
p     = xout(:,11) * 180/pi;
q     = xout(:,12) * 180/pi;
r     = xout(:,13) * 180/pi; 

%% plots
figure(1)
subplot(211)
plot(y,x,'linewidth',2)
grid,axis('equal'),xlabel('East'),ylabel('North'),title('Semisub xy-plot (m)')
subplot(212)
plot(t,x,t,y,'linewidth',2),title('Semisub positions (m)'),
legend('x position','y position'),grid

figure(2)
subplot(311),plot(t,phi,'linewidth',2),xlabel('time (s)')
title('roll angle \phi (deg)'),grid
subplot(312),plot(t,theta,'linewidth',2),xlabel('time (s)')
title('pitch angle \theta (deg)'),grid
subplot(313),plot(t,psi,'linewidth',2),xlabel('time (s)')
title('yaw angle \psi (deg)'),grid


