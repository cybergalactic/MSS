function [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t_k)
% insSignal is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an INS signal generator for testing of Inertial
% Navigation Systems (INS).
%
% Syntax:
%   [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t_k)
%
% Inputs:
%   x[k]: A 15-element signal vector structured as follows:
%         p^n - NED (North-East-Down) position vector, 3 elements.
%         v^n - NED velocity vector, 3 elements.
%         b_acc^b - Acceleration bias in the BODY frame, 3 elements.
%         Theta - Euler angles (phi, theta, psi), 3 elements.
%         b_ars^n - Attitude Rate Bias (ARS) in the BODY frame, 3 elements.
%
%   mu   : Latitude in radians.
%   h    : Sampling time in seconds.
%   t_k  : Current time t[k] in seconds.
%
% Outputs:
%   x[k+1]     : Propagated signal vector, 15 elements
%   f_imu[k+1] : IMU specific force vector expressed in the BODY frame, 
%                units in meters/second^2.
%   w_imu[k+1] : IMU attitude rate vector expressed in the BODY frame, 
%                units in radians/second.
%   m_imu[k+1] : IMU magnetic field vector expressed in the BODY frame, 
%                units in nanoTesla (nT).
%   m_ref      : Reference magnetic field vector expressed in the NED frame, 
%                specific to Trondheim, Norway, units in nT.
%
% Author:
%   Thor I. Fossen
% Date:
%   2020-04-28
% Revisions:
%   2024-08-23 : Use RK4 for state propagation, update signal generators

m_ref = [13559 921 50209]'; % Reference magnetic field, Trondheim, Norway
g_n = [0 0 gravity(mu)]'; % NED gravity vector for specified lattitude mu

% True specific force: f[k] = a[k] - g
f_true = @(t, g_n, x) [
    0.1 * sin(0.1 * t);  
    0.1 * cos(0.1 * t);
    0.05 * sin(0.05 * t)] - Rzyx(x(10),x(11),x(12))' * g_n;

% True angular rate: w[k]
w_true = @(t) [
    0.02 * cos(0.2 * t);
   -0.02 * sin(0.1 * t);
    0.01 * sin(0.1 * t) ] ;

% Define vehicle dynamics as an inline function
dynamics = @(x, t) [
    x(4:6); % Linear velocity 
    Rzyx(x(10),x(11),x(12)) * f_true(t,g_n,x) + g_n; % Linear acceleration 
    zeros(3,1); % Acceleration bias (static in this model)
    Tzyx(x(10),x(11)) * w_true(t); % Angular rates 
    zeros(3,1)]; % Angular rate bias (static in this model)

% Propagate states: x[k+1]
x = rk4(dynamics, h, x, t_k);

% IMU inertial measurements at time t_k+1 with noise addition
w1 = 0.01 * randn(3,1); % Gaussian white noise for specific force
w2 = 0.01 * randn(3,1); % Gaussian white noise for angular rate
b_acc = x(7:9); % ACC bias: b_acc[k+1]
b_ars = x(13:15); % ARS bias: b_ars[k+1]
f_imu = f_true(t_k+h, g_n, x) + b_acc + w1; % IMU specific force: f_imu[k+1]
w_imu = w_true(t_k+h) + b_ars + w2; % IMU angular rate: w_imu[k+1]

% Calculation of magnetic field readings
w3 = 0.01 * randn(3,1); % Gaussian white noise for magnetic field
m_imu = Rzyx(x(10),x(11),x(12))' * m_ref + w3; % IMU magnetic field: m_imu[k+1]

end
