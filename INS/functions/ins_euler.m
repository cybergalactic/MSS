function [x_ins, P_prd] = ins_euler(...
    x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, psi, y_pos, y_vel)
% ins_euler is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter
% (ESKF) specifically for Inertial Navigation Systems (INS) that are aided
% by compass and positional data. Attitude is parametrized using the 3-parameter
% Euler angle representation, which is singular for theta = +- 90 deg. To avoid 
% the singularity, use ins_mekf or ins_mekf_psi. instead.
%
% Usage scenarios are detailed in examples SIMaidedINSeuler and exINS_Euler,
% demonstrating the implementation of the Kalman filter loop using the
% corrector-predictor representation:
%
%   - With new slow position measurements:
%       [x_ins,P_prd] = ins_euler(...
%           x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos)
%       [x_ins,P_prd] = ins_euler(...
%           x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos, y_vel)
%
%   - Without new positions measurements (no aiding):
%       [x_ins,P_prd] = ins_euler(...
%          x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi)
%
% This function models the INS errors in a 15-dimensional state space,
% including position, velocity, biases, and attitude errors:
%
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k]
%
% Inputs:
%   x_ins[k] : INS state vector at step k, includes position, velocity,
%              accelerometer biases, attitude (Euler angles), and gyro biases.
%   P_prd[k] : 15x15 covariance matrix of the prediction step.
%   mu       : Latitude in radians, used to calculate Earth's gravity vector.
%   h        : Sampling time in seconds.
%   Qd, Rd   : Process and measurement noise covariance matrices for the
%              Kalman filter.
%   f_imu[k] : Specific force measurements from the IMU.
%   w_imu[k] : Angular rate measurements from the IMU.
%   y_psi[k] : Fast compass measurement (yaw angle).
%   y_pos[k] : Slow position measurements aids the filter.
%   y_vel[k] : (Optionally) Slow velocity measurements aids the filter.
%
% Outputs:
%   x_ins[k+1] - Updated INS state vector after propagation.
%   P_prd[k+1] - Updated prediction covariance matrix after propagation.
%
% References:
%   T. I. Fossen (2021). "Handbook of Marine Craft Hydrodynamics and Motion
%   Control," 2nd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2021-01-14
% Revisions:
%   2021-12-21: Improved numerical accuracy by replacing Euler's method
%               with exact discretization in the INS PVA propagation.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method
%               with RK4 in the INS attitude dynamics.
%   2024-08-31: Using invQR.m instead of inv.m

% Bias time constants (user specified)
T_acc = 1000;
T_ars = 500;

%% ESKF states and matrices
p_ins = x_ins(1:3);          % INS states
v_ins = x_ins(4:6);
b_acc_ins = x_ins(7:9);
theta_ins = x_ins(10:12);
b_ars_ins = x_ins(13:15);

% Gravity vector
g = gravity(mu);            % WGS-84 gravity model
g_n = [0 0 g]';

% Constants
O3 = zeros(3,3);
I3 = eye(3);

% Rotation matrix
R = Rzyx(theta_ins(1), theta_ins(2), theta_ins(3));

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Normalized specific force 
v01 = [0 0 -1]';          % NED gravity reference vector, f_z = [0 0 -g]' at rest
v1 = f_ins / norm(f_ins); % BODY specific force measurement

% Discrete-time ESKF matrices
A = [ O3 I3  O3            O3              O3
      O3 O3 -R            -R*Smtrx(f_ins)  O3
      O3 O3 -(1/T_acc)*I3  O3              O3
      O3 O3  O3           -Smtrx(w_ins)   -I3
      O3 O3  O3            O3             -(1/T_ars)*I3 ];

% Ad = eye(15) + h * A + 0.5 * (h * A)^2 + ...
Ad = expm_taylor(A * h); 

if (nargin == 10)
    Cd = [ I3 O3 O3 O3 O3                % NED positions
           O3 O3 O3 Smtrx(R'*v01) O3     % Gravity
           zeros(1,11) 1 zeros(1,3) ];   % Compass
else
    Cd = [ I3 O3 O3 O3 O3                % NED positions
           O3 I3 O3 O3 O3                % NED velocities
           O3 O3 O3 Smtrx(R'*v01) O3     % Gravity
           zeros(1,11) 1 zeros(1,3)  ];  % Compass
end

Ed = h *[ O3 O3    O3 O3
         -R  O3    O3 O3
          O3 I3    O3 O3
          O3 O3   -I3 O3
          O3 O3    O3 I3  ];

%% Kalman filter algorithm
if (nargin == 9)             % No aiding

    P_hat = P_prd;

else                         % INS aiding

    % ESKF gain: K[k]
    K = P_prd * Cd' * invQR(Cd * P_prd * Cd' + Rd);
    IKC = eye(15) - K * Cd;

    % Estimation error: eps[k]
    eps_pos = y_pos - p_ins;
    eps_g   = v1 - R' * v01;
    eps_psi = ssa(psi - theta_ins(3));

    if (nargin == 10)
        eps = [eps_pos; eps_g; eps_psi];
    else
        eps_vel = y_vel - v_ins;
        eps = [eps_pos; eps_vel; eps_g; eps_psi];
    end

    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * eps;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';

    % INS reset: x_ins[k]
    p_ins = p_ins + delta_x_hat(1:3);	         % Reset INS position
    v_ins = v_ins + delta_x_hat(4:6);			 % Reset INS velocity
    b_acc_ins = b_acc_ins + delta_x_hat(7:9);    % Reset ACC bias
    theta_ins = theta_ins + delta_x_hat(10:12);  % Reset INS attitude
    b_ars_ins = b_ars_ins + delta_x_hat(13:15);  % Reset ARS bias

end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: x_ins[k+1]
a_ins = R * f_ins + g_n;                        % Linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;      % Exact discretization
v_ins = v_ins + h * a_ins;                      % Exact discretization
theta_ins = rk4(@attitudeDynamics , h, theta_ins, w_ins); % RK4

% Euler's method, alternative to rk4
% theta_ins = theta_ins + h * Tzyx(theta_ins(1), theta_ins(2)) * w_ins;         

x_ins = [p_ins; v_ins; b_acc_ins; theta_ins; b_ars_ins];

%% Attitude dynamics function
function thetadot = attitudeDynamics(theta, w)
    % attitudeDynamics computes the time derivative of the Euler angles.
    thetadot = Tzyx(theta(1), theta(2)) * w;
end

end