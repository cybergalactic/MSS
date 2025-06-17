function [x, f_imu, w_imu, m_imu] = insSignal(x, h, t_k, mu, m_ref, signalNo)
% insSignal is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an INS signal generator for testing of Inertial
% Navigation Systems (INS).
%
% Syntax:
%   [x, f_imu, w_imu, m_imu] = insSignal(x, h, t_k, mu, m_ref)
%   [x, f_imu, w_imu, m_imu] = insSignal(x, h, t_k, mu, m_ref, signalNo)
%
% Inputs:
%   x[k]  : A 15-element signal vector structured as follows:
%           p^n - NED (North-East-Down) position vector, 3 elements.
%           v^n - NED velocity vector, 3 elements.
%           b_acc^b - Acceleration bias in the BODY frame, 3 elements.
%           Theta - Euler angles (phi, theta, psi), 3 elements.
%           b_ars^n - Attitude Rate Bias (ARS) in the BODY frame, 3 elements.
%   h     : Sampling time in seconds.
%   t_k   : Current time t_k in seconds.
%   mu    : Latitude in radians, see magnetifField.m.
%   m_ref : Reference magnetic field vector expressed in NED,
%           see magnetifField.m.
%   signalNo : Integer for choosing test signal, default 1 if argument is omitted
%
% Outputs:
%   x[k+1]     : Propagated signal vector, 15 elements
%   f_imu[k+1] : IMU specific force vector expressed in the BODY frame, 
%                units in meters/second^2.
%   w_imu[k+1] : IMU attitude rate vector expressed in the BODY frame, 
%                units in radians/second.
%   m_imu[k+1] : IMU magnetic field vector expressed in the BODY frame, 
%                units in nanoTesla (nT).
%
% Author:
%   Thor I. Fossen
% Date:
%   2020-04-28
% Revisions:
%   2024-08-23 : New improved version of the signal generator.

if nargin == 5
    signalNo = 1;
end

g_n = [0 0 gravity(mu)]'; % NED gravity vector for lattitude mu

switch signalNo
    case 1
        % True specific force: f[k] = a[k] - g
        f_true = @(t, g_n, x) [
            0.1 * sin(0.1 * t);
            0.01 * cos(0.1 * t);
            0.05 * sin(0.05 * t)
            ] - Rzyx(x(10), x(11), x(12))' * g_n;
        
        % True angular rate: w[k]
        w_true = @(t) [
            0.03 * cos(0.2 * t);
            -0.02 * sin(0.1 * t);
            0.01 * sin(0.05 * t) 
            ];
    case 2
        % True specific force: f[k] = a[k] - g
        f_true = @(t, g_n, x) [
            0.2 * sin(0.07 * t + 0.5) + 0.05 * sin(0.2 * t);
            0.1 * cos(0.1 * t) + 0.02 * cos(0.3 * t + 0.3);
            0.15 * sin(0.05 * t + 1.0) + 0.03 * cos(0.2 * t)
            ] - Rzyx(x(10), x(11), x(12))' * g_n;

        % True angular rate: w[k]
        w_true = @(t) [
            0.05 * sin(0.3 * t) + 0.01 * cos(0.1 * t + 0.7);
            0.04 * cos(0.2 * t + 1.2);
            0.03 * sin(0.15 * t + 0.4)
            ];
    otherwise
        disp('Use integer 1 or 2 for test signal');

end

% Define vehicle dynamics as an inline function
dynamics = @(x, t) [
    x(4:6); % Linear velocity 
    Rzyx(x(10),x(11),x(12)) * f_true(t,g_n,x) + g_n; % Linear acceleration 
    zeros(3,1); % Acceleration bias (static in this model)
    Tzyx(x(10),x(11)) * w_true(t); % Angular rates 
    zeros(3,1) ]; % Angular rate bias (static in this model)

% Propagate states: x[k+1]
x = rk4(dynamics, h, x, t_k);

% IMU inertial measurements at time t_k+1 with noise addition
w1 = 0.005 * randn(3,1); % Gaussian white noise for specific force
w2 = 0.005 * randn(3,1); % Gaussian white noise for angular rate
b_acc = x(7:9); % ACC bias: b_acc[k+1]
b_ars = x(13:15); % ARS bias: b_ars[k+1]
f_imu = f_true(t_k+h, g_n, x) + b_acc + w1; % IMU specific force: f_imu[k+1]
w_imu = w_true(t_k+h) + b_ars + w2; % IMU angular rate: w_imu[k+1]

% Calculation of magnetic field readings
w3 = 0.005 * randn(3,1); % Gaussian white noise for magnetic field
m_imu = Rzyx(x(10),x(11),x(12))' * m_ref + w3; % IMU magnetic field: m_imu[k+1]

end
