% SIMrig    User editable script for simulation of a semi-submersible with
% two pontoons and four legs. The data is generated using rig.m, which
% calls data_rig.m. 
%
%   m =  27 162 500 (mass)
%   H_p = 13.0 (pontoon height)
%   B_p = 12.0 (pontoon width)
%   L_p = 84.6 (pontoon length)
%   r_CG = [-0.5 0 -19.5]' (center of gravity in COO)
%
% Calls:       euler2.m, Rzyx.m, eulerang.m
%
% Author:      Thor I. Fossen
% Date:        2019-03-12
% Revisions: 

clear all

f = 20;              % sampling frequency
h = 1/f;             % sampling time
N = round(1000*f);   % number of samples

w_n = 0.02;          % Bandwidth
Kp = w_n^2;          % Proportional gain 
Kd = 2 *w_n^2;       % Derivative  gain 
Ki = 0.1 * Kp * w_n; % Integral  gain 

% Regular wave data
beta_wave = 60;  % wave direction 
A_wave = 2000000; % wave gain 
X_wave = A_wave * cos(beta_wave * pi/180);
Y_wave = A_wave * sin(beta_wave * pi/180);

% initial values
eta = zeros(6,1);   % generalized position
nu  = zeros(6,1);   % generalized velcoity
z_3 = zeros(3,1);   % integral state

% Semi-sub model
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

M_33 = [M(1,1) M(1,2) M(1,6)          % 3-DOF inertia mattrix
        M(2,1) M(2,2) M(2,6)
        M(6,1) M(6,2) M(6,6) ];

    
% *************** MAIN SIMULATION LOOP ************************
for i=1:N+1
    t = (i-1)*h;                   % simulation time in seconds

    % kinematics
    phi = eta(4); theta = eta(5); psi = eta(6);
    J = eulerang(phi,theta,psi); 
    R = Rzyx(phi,theta,psi); 
    
    % 3-DOF nonlinear DP control system
    eta_3 = [eta(1) eta(2) eta(6)]' ;
    nu_3  = [nu(1) nu(2) nu(6)]';   
    tau_3 = -M_33 * (Kp * R' * eta_3 - Kd * nu_3 + Ki * R' * z_3); 
    tau = [tau_3(1) tau_3(2) 0 0 0 tau_3(3)]';
    
    % rig model
    t_wave = zeros(6,1);
    t_wave(1) = X_wave * sin(0.2 * t);
    t_wave(2) = Y_wave * sin(0.2 * t);   
    etadot = J * nu;
    nudot = Minv * (tau - D * nu + t_wave);
    
    % store data for presentation
    xout(i,:) = [t, eta', nu'];
    
    % numerical integration
    eta = euler2(etadot, eta, h);             % Euler integration
    nu  = euler2(nudot, nu, h);     
    z_3 = euler2(eta_3, z_3, h);

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

% plots
figure(1)
plot(y,x)
grid,axis('equal'),xlabel('East'),ylabel('North'),title('Semi-sub position (m)')

figure(2)
subplot(311),plot(t,phi),xlabel('time (s)'),title('roll angle \phi (deg)'),grid
subplot(312),plot(t,theta),xlabel('time (s)'),title('pitch angle \theta (deg)'),grid
subplot(313),plot(t,psi),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid


