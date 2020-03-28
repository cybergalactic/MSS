function [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t_k)
% INS signal generator
%
%    [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t_k)
%
% Inputs: 
%    x, signal vector [p^n; v^n; b_acc^b; Theta; b_ars^n], dim 15 
%       
%       p^n, NED position vector, dim 3 
%       v^n, NED velocity vector, dim 3 
%       b_acc^n, BODY acceleration bias vector, dim 3 
%       Theta, vector of Euler angles phi, theta and psi, dim 3 
%       b_ars^n, BODY attitude rate bias vector, dim 3 
%
%    mu, lattitude [rad]
%    h, sampling time [s]
%    t_k, present time [s] 
%
% Outputs: 
%    x, propgated signal vector, dim 15
%    f_imu, IMU specific force vector expressed in BODY [m/s^2]
%    w_imu, IMU attitude rate vector expressed in BODY [rad/s]
%    m_imu, IMU magentic field vector expressed in BODY [nT]
%    m_ref, Magentic field vector expressed in NED for Trondheim, Norway
%
% Author:    Thor I. Fossen
% Date:      28 March 2020
% Revisions: 

g_n = [ 0 0 gravity(mu) ]';             % NED gravity vector
g_b = Rzyx(x(10),x(11),x(12))' * g_n;   % BODY gravity vector
    
f_true = [0.1  * sin(0.1*t_k)           % true specific force: f = a - g
          0.1  * cos(0.1*t_k)
         0.05 * sin(0.05*t_k) ] - g_b;

w_true = [ 0.01 * cos(0.2*t_k)              % true angular rate
          -0.02 * sin(0.1*t_k)
           0.01 * sin(0.1*t_k) ];

x_dot = [ x(4:6)                            % true states
    Rzyx(x(10),x(11),x(12)) * f_true + g_n
    zeros(3,1)
    Tzyx(x(10),x(11)) * w_true
    zeros(3,1) ];

% IMU measurements
w1 = 0.01 * randn(3,1);                      % white noise
w2 = 0.01 * randn(3,1);
b_acc = x(7:9);
b_ars = x(13:15);
f_imu = f_true + b_acc + w1;
w_imu = w_true + b_ars + w2;

% NOA magentic field calculator 
% https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
w3 = 0.01 * randn(3,1);       % white noise
m_ref = [13559 921 50209]';   % [nT] - Location: Trondheim, Norway
m_imu = Rzyx(x(10),x(11),x(12))' * m_ref + w3;

% Propagate signal vector
x = x + h * x_dot;

end   

