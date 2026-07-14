function [x, f_imu, w_imu, m_imu] = insSignal(x, h, t_k, mu, m_ref, signalNo)
% insSignal is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an INS signal generator for testing of Inertial
% Navigation Systems (INS).
%
% Syntax:
%   [x, f_imu, w_imu, m_imu] = insSignal(x, h, t_k, mu, m_ref)
%   [x, f_imu, w_imu, m_imu] = insSignal( ...
%       x, h, t_k, mu, m_ref, signalNo)
%
% Inputs:
%   x[k]     : A 15-element signal vector structured as follows:
%              p^n       - NED position vector, 3 elements
%              v^n       - NED velocity vector, 3 elements
%              b_acc^b   - Accelerometer bias in BODY, 3 elements
%              Theta     - Euler angles (phi, theta, psi), 3 elements
%              b_ars^b   - ARS bias in BODY, 3 elements
%   h        : Sampling time [s]
%   t_k      : Current time [s]
%   mu       : Latitude [rad], see magneticField.m
%   m_ref    : Reference magnetic field vector expressed in NED
%   signalNo : Test-signal selector:
%              1 - Constant IMU biases
%              2 - Slowly varying Gauss--Markov IMU biases
%
% Outputs:
%   x[k+1]     : Propagated signal vector, 15 elements
%   f_imu[k+1] : IMU specific force in BODY [m/s^2]
%   w_imu[k+1] : IMU angular rate in BODY [rad/s]
%   m_imu[k+1] : IMU magnetic field in BODY [nT]
%
% Author:
%   Thor I. Fossen
% Date:
%   2020-04-28
% Revisions:
%   2024-08-23: New improved version of the signal generator.
%   2026-07-14: Added slowly varying Gauss--Markov biases for test signal 2.

% Default test signal #1 (constant biases)
if nargin == 5
    signalNo = 1;
end

% Bias-model parameters used for test signal 2
T_acc = 300;   % Accelerometer-bias correlation time [s]
T_ars = 300;   % ARS-bias correlation time [s]

sigma_b_acc = 20e-6 * 9.81;     % Stationary standard deviation [m/s^2]
sigma_b_ars = deg2rad(8/3600);  % Stationary standard deviation [rad/s]

% NED gravity vector
g_n = [0 0 gravity(mu)]';       

switch signalNo

    case 1

        % True specific force: f = a - g
        f_true = @(t, g_n, x) [
            0.1  * sin(0.1  * t);
            0.01 * cos(0.1  * t);
            0.05 * sin(0.05 * t)
            ] - Rzyx(x(10), x(11), x(12))' * g_n;

        % True angular rate
        w_true = @(t) [
             0.03 * cos(0.2  * t);
            -0.02 * sin(0.1  * t);
             0.01 * sin(0.05 * t)
            ];

    case 2

        % True specific force: f = a - g
        f_true = @(t, g_n, x) [
            0.2  * sin(0.07 * t + 0.5) ...
                + 0.05 * sin(0.2 * t);
            0.1  * cos(0.1 * t) ...
                + 0.02 * cos(0.3 * t + 0.3);
            0.15 * sin(0.05 * t + 1.0) ...
                + 0.03 * cos(0.2 * t)
            ] - Rzyx(x(10), x(11), x(12))' * g_n;

        % True angular rate
        w_true = @(t) [
            0.05 * sin(0.3  * t) ...
                + 0.01 * cos(0.1 * t + 0.7);
            0.04 * cos(0.2  * t + 1.2);
            0.03 * sin(0.15 * t + 0.4)
            ];

    otherwise
        error('signalNo must be 1 or 2.');

end

% Vehicle dynamics. Biases are propagated separately below.
dynamics = @(x, t) [
    x(4:6);
    Rzyx(x(10), x(11), x(12)) * f_true(t, g_n, x) + g_n;
    zeros(3,1);
    Tzyx(x(10), x(11)) * w_true(t);
    zeros(3,1)
    ];

% Propagate position, velocity, and attitude
x = rk4(dynamics, h, x, t_k);

% Bias propagation: For signalNo == 1, x(7:9) and x(13:15) remain constant
if signalNo == 2

    % Exact discrete-time first-order Gauss--Markov model
    alpha_acc = exp(-h / T_acc);
    alpha_ars = exp(-h / T_ars);

    x(7:9) = alpha_acc * x(7:9) ...
        + sigma_b_acc * sqrt(1 - alpha_acc^2) * randn(3,1);

    x(13:15) = alpha_ars * x(13:15) ...
        + sigma_b_ars * sqrt(1 - alpha_ars^2) * randn(3,1);
end

% IMU measurements at time t_k + h
w1 = 0.005 * randn(3,1);   % Specific-force measurement noise
w2 = 0.005 * randn(3,1);   % Angular-rate measurement noise

b_acc = x(7:9);
b_ars = x(13:15);

f_imu = f_true(t_k + h, g_n, x) + b_acc + w1;
w_imu = w_true(t_k + h) + b_ars + w2;

% Magnetometer measurement
w3 = 0.005 * randn(3,1);
m_imu = Rzyx(x(10), x(11), x(12))' * m_ref + w3;

end