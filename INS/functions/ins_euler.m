function [x_ins, P_prd] = ins_euler( ...
    x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, imu_meas, y_psi, y_pos, y_vel)
% ins_euler is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter
% (ESKF) for an Inertial Navigation System (INS) aided by compass and
% position measurements. Velocity aiding is optional.
%
% Attitude is parameterized using the three-parameter Euler-angle
% representation, which is singular for theta = +/- 90 deg. To avoid this
% singularity, use quaternions; see ins_mekf or ins_mekf_psi.
%
% Usage scenarios are detailed in SIMaidedINSeuler demonstrating the 
% corrector-predictor implementation:
%
%   - With new position measurements:
%       [x_ins,P_prd] = ins_euler( ...
%           x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, ...
%           imu_meas, y_psi, y_pos)
%
%   - With new position and velocity measurements:
%       [x_ins,P_prd] = ins_euler( ...
%           x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, ...
%           imu_meas, y_psi, y_pos, y_vel)
%
%   - Without new position or velocity measurements:
%       [x_ins,P_prd] = ins_euler( ...
%           x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, ...
%           imu_meas, y_psi)
%
% The 15-dimensional error-state vector is
%
%   delta_x = [delta_p; delta_v; delta_b_acc;
%              delta_theta; delta_b_ars]
%
% with discrete-time model
%
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k].
%
% Inputs:
%   x_ins[k]    : 15-element INS state vector containing position, velocity,
%                 accelerometer bias, Euler angles, and ARS bias.
%   P_prd[k]    : 15x15 predicted error covariance matrix.
%   mu          : Latitude in radians, used to calculate the gravity vector.
%   h           : Sampling time in seconds.
%   Qd          : Process-noise covariance matrix.
%   Rd          : Measurement-noise covariance matrix.
%   T_acc       : Acceleration bias time constant in seconds.
%   T_ars       : Angular rate bias time constant in seconds.
%   imu_meas[k] : Six-element vector containing [f_imu; w_imu].
%                 [fx, fy, fz, wx, wy, wz]. More specific, f_imu[k] is a 3x1 
%                  vector representing the IMU specific force measurements and 
%                  w_imu[k] is a 3x1 vector representing the IMU angular rate 
%                  measurements. The IMU axes are assumed to be oriented 
%                  forward-starboard-down.
%   y_psi[k]    : Compass heading measurement.
%   y_pos[k]    : Position measurement expressed in NED.
%   y_vel[k]    : Optional velocity measurement expressed in NED.
%
% Outputs:
%   x_ins[k+1]  : Updated and propagated INS state vector.
%   P_prd[k+1]  : Predicted error covariance matrix.
%
% References:
%   T. I. Fossen (2027). "Handbook of Marine Craft Hydrodynamics and Motion
%   Control," 3rd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2021-01-14
% Revisions:
%   2021-12-21: Improved numerical accuracy by replacing Euler's method
%               with constant-acceleration INS PVA propagation.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method
%               with RK4 in the INS attitude dynamics.
%   2026-07-15: Added consistent additive Euler-angle error dynamics and
%               reference-vector measurement Jacobian.

% ==============================================================================
% INS states and constants
% ==============================================================================
p_ins     = x_ins(1:3);
v_ins     = x_ins(4:6);
b_acc_ins = x_ins(7:9);
theta_ins = x_ins(10:12);
b_ars_ins = x_ins(13:15);

% WGS-84 gravity vector expressed in NED
g_n = [0 0 gravity(mu)]'; 

% Matrix constants
O3 = zeros(3,3);
I3 = eye(3);

% ==============================================================================
% Quantities evaluated at the predicted INS state
% ==============================================================================

% Euler-angle kinematic matrices
R = Rzyx(theta_ins(1), theta_ins(2), theta_ins(3));
T = Tzyx(theta_ins(1), theta_ins(2));

% Bias-compensated IMU measurements
imu_meas = imu_meas(:);
f_imu = imu_meas(1:3);
w_imu = imu_meas(4:6);

f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Normalized specific-force reference-vector measurement
v01 = [0 0 -1]';              % NED gravity reference vector
v1  = f_ins / norm(f_ins);    % BODY specific-force measurement

% ==============================================================================
% Discrete-time measurement matrix
% ==============================================================================
if nargin == 11

    % Position, gravity reference vector, and compass
    Cd = [ ...
        I3           O3 O3 O3                    O3;           % NED position
        O3           O3 O3 Smtrx(R' * v01) / T   O3;           % Gravity vector
        zeros(1,11)  1                           zeros(1,3) ]; % Compass heading

elseif nargin == 12

    % Position, velocity, gravity reference vector, and compass
    Cd = [ ...
        I3            O3 O3 O3     O3;                % NED position
        O3            I3 O3 O3     O3;                % NED velocity
        O3            O3 O3 Smtrx(R' * v01) / T   O3; % Gravity vector
        zeros(1,11)   1  zeros(1,3) ];                % Compass heading

end

% ==============================================================================
% Kalman filter corrector
% ==============================================================================
if nargin == 10

    % No new position or velocity aiding measurements
    P_hat = P_prd;

else

    % ESKF gain: K[k]
    K = P_prd * Cd' / (Cd * P_prd * Cd' + Rd);
    IKC = eye(15) - K * Cd;

    % Measurement innovations
    eps_pos = y_pos - p_ins;
    eps_g   = v1 - R' * v01;
    eps_psi = ssa(y_psi - theta_ins(3));

    if nargin == 11
        eps = [eps_pos; eps_g; eps_psi];
    else
        eps_vel = y_vel - v_ins;
        eps = [eps_pos; eps_vel; eps_g; eps_psi];
    end

    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * eps;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';

    % Feedback reset of the INS states
    p_ins = p_ins + delta_x_hat(1:3);
    v_ins = v_ins + delta_x_hat(4:6);
    b_acc_ins = b_acc_ins + delta_x_hat(7:9);
    theta_ins = theta_ins + delta_x_hat(10:12);
    b_ars_ins = b_ars_ins + delta_x_hat(13:15);

    % Recompute quantities affected by the feedback reset
    R = Rzyx(theta_ins(1), theta_ins(2), theta_ins(3));
    T = Tzyx(theta_ins(1), theta_ins(2));

    f_ins = f_imu - b_acc_ins;
    w_ins = w_imu - b_ars_ins;

end

% ==============================================================================
% ESKF covariance predictor
% ==============================================================================
Atheta = eulerRateJacobian(theta_ins,w_ins); % Jacobian of the Euler-angle kinematics

% Continuous-time additive Euler-angle error dynamics
A = [ ...
    O3 I3  O3               O3                         O3;
    O3 O3 -R               -R * Smtrx(f_ins) / T       O3;
    O3 O3 -(1/T_acc) * I3   O3                         O3;
    O3 O3  O3               Atheta                     -T;
    O3 O3  O3               O3                        -(1/T_ars) * I3 ];

% Discrete-time state-transition matrix
Ad = expm_taylor(A * h);

% Discrete-time process-noise input matrix
Ed = h * [ ...
     O3 O3  O3 O3;
     -R O3  O3 O3;
     O3 I3  O3 O3;
     O3 O3  -T O3;
     O3 O3  O3 I3 ];

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% ==============================================================================
% INS propagation
% ==============================================================================
a_ins = R * f_ins + g_n; % Navigation-frame acceleration
p_ins = p_ins + h * v_ins + 0.5 * h^2 * a_ins;
v_ins = v_ins + h * a_ins;
theta_ins = rk4(@attitudeDynamics, h, theta_ins, w_ins);

% Euler's method, alternative to RK4:
% theta_ins = theta_ins + h * Tzyx(theta_ins(1), theta_ins(2)) * w_ins;

% Updated INS state vector
x_ins = [ ...
    p_ins;
    v_ins;
    b_acc_ins;
    theta_ins;
    b_ars_ins ];


% ==============================================================================
%% FUNCTION: Attitude dynamics 
% ==============================================================================
function theta_dot = attitudeDynamics(theta, w)
    theta_dot = Tzyx(theta(1), theta(2)) * w; % Time derivative of the Euler angles
end

% ==============================================================================
%% FUNCTION: Jacobian of Tzyx(theta) * w with respect to the Euler angles
% ==============================================================================
function Atheta = eulerRateJacobian(theta, w)
    phi = theta(1);
    th  = theta(2);
    q = w(2);
    r = w(3);

    sphi = sin(phi); 
    cphi = cos(phi); 
    tth = tan(th); 
    cth = cos(th); 
    secth = 1 / cth;

    % Partial derivative with respect to phi
    d_dphi = [ ...
        tth   * (cphi*q - sphi*r);
                -sphi*q - cphi*r;
        secth * (cphi*q - sphi*r) ];

    % Partial derivative with respect to theta
    d_dtheta = [ ...
        secth^2   * (sphi*q + cphi*r);
        0;
        secth*tth * (sphi*q + cphi*r) ];

    % Tzyx does not depend explicitly on psi
    Atheta = [ ...
        d_dphi, ...
        d_dtheta, ...
        zeros(3,1) ];
end

end