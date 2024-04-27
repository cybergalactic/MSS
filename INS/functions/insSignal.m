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
%   mu  : Latitude in radians.
%   h   : Sampling time in seconds.
%   t_k : Current time in seconds.
%
% Outputs:
%   x[k1] : Signal vector, 15 elements.
%   f_imu : IMU specific force vector expressed in the BODY frame, units in meters/second^2.
%   w_imu : IMU attitude rate vector expressed in the BODY frame, units in radians/second.
%   m_imu : IMU magnetic field vector expressed in the BODY frame, units in nanoTesla (nT).
%   m_ref : Reference magnetic field vector expressed in the NED frame, specific to Trondheim, Norway, units in nT.
%
% Author:
%   Thor I. Fossen
% Date:
%   2020-04-28
% Revisions:
%   None

g_n = [0 0 gravity(mu)]';                 % NED gravity vector
g_b = Rzyx(x(10),x(11),x(12))' * g_n;     % Transforming the gravity vector 
                                          % from NED to BODY

f_true = [0.1 * sin(0.1*t_k)              % True specific force: f = a - g
          0.1 * cos(0.1*t_k)
          0.05 * sin(0.05*t_k)] - g_b;

w_true = [0.02 * cos(0.2*t_k)             % True angular rate
          -0.02 * sin(0.1*t_k)
           0.01 * sin(0.1*t_k)];

x_dot = [x(4:6)                           % Derivative of true states
         Rzyx(x(10),x(11),x(12)) * f_true + g_n
         zeros(3,1)
         Tzyx(x(10),x(11)) * w_true
         zeros(3,1)];

% Simulation of IMU measurements with noise addition
w1 = 0.01 * randn(3,1);         % Gaussian white noise for specific force
w2 = 0.01 * randn(3,1);         % Gaussian white noise for angular rate
b_acc = x(7:9);
b_ars = x(13:15);
f_imu = f_true + b_acc + w1;   % IMU specific force output
w_imu = w_true + b_ars + w2;   % IMU angular rate output

% Calculation of magnetic field readings
w3 = 0.01 * randn(3,1);        % Gaussian white noise for magnetic field
m_ref = [13559 921 50209]';    % Reference magnetic field, Trondheim, Norway
m_imu = Rzyx(x(10),x(11),x(12))' * m_ref + w3; % IMU magnetic field output

% Propagation of the signal vector
x = x + h * x_dot;

end
