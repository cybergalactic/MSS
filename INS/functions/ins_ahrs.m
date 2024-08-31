function [x_ins, P_prd] = ins_ahrs( ...
    x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_ahrs, y_pos, y_vel)
% ins_ahrs is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter 
% (ESKF) specifically for Inertial Navigation Systems (INS) that are 
% augmented by an attitude heading reference systems (AHRS) and aided by 
% positional data. Attitude is parametrized using the 3-parameter Euler 
% angle representation, which is singular for theta = +- 90 deg.
%
% Usage scenarios are detailed in examples SIMaidedINSeuler and ExINS_AHRS, 
% demonstrating the implementation of the Kalman filter loop using the 
% corrector-predictor representation:
%
%   - With new slow position measurements:
%       [x_ins,P_prd] = ins_ahrs(...
%           x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_ahrs, y_pos)
%       [x_ins,P_prd] = ins_ahrs(...
%           x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_ahrs, y_pos, y_vel)
%
%   - Without new position measurements (no aiding):
%       [x_ins,P_prd] = ins_ahrs(x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs)
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
%   y_ahrs[k]: Attitude measurements (roll, pitch, yaw) from the AHRS.
%   y_pos[k] : Slow position measurements aids the filter.
%   y_vel[k] : (Optionally) Slow velocity measurements aids the filter.
%
% Outputs:
%   x_ins[k+1] : Updated INS state vector after propagation.
%   P_prd[k+1] : Updated prediction covariance matrix after propagation.
%
% References:
%   T. I. Fossen (2021). "Handbook of Marine Craft Hydrodynamics and Motion 
%   Control," 2nd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2020-03-21
% Revisions: 
%   2021-12-21: Improved numerical accuracy by replacing Euler's method
%               with exact discretization in the INS PVA propagation.
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
g_n = [0 0 gravity(mu)]';    % WGS-84 gravity model

% Constants 
O3 = zeros(3,3);
I3 = eye(3);

% Transformation matrices
R = Rzyx(y_ahrs(1), y_ahrs(2), y_ahrs(3));
T = Tzyx(y_ahrs(1), y_ahrs(2)); 

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Discrte-time ESKF matrices
A = [ O3 I3  O3           O3  O3
      O3 O3 -R            O3  O3
      O3 O3 -(1/T_acc)*I3 O3  O3
      O3 O3  O3           O3 -T
      O3 O3  O3           O3 -(1/T_ars)*I3 ];
   
% Ad = eye(15) + h * A + 0.5 * (h * A)^2 + ...
Ad = expm_taylor(A * h); 

if (nargin == 10)
    Cd = [ I3 O3 O3 O3 O3        % NED positions (x, y, z)
           O3 O3 O3 I3 O3];      % Euler angles (phi, theta, psi)
else
    Cd = [ I3 O3 O3 O3 O3        % NED positions (x, y, z)
           O3 I3 O3 O3 O3        % NED velocities       
           O3 O3 O3 I3 O3];      % Euler angles (phi, theta, psi)
end
       
Ed = h *[  O3 O3    O3 O3
          -R  O3    O3 O3
           O3 I3    O3 O3
           O3 O3   -T  O3
           O3 O3    O3 I3  ];

%% Kalman filter algorithm       
if (nargin == 9)             % No aiding
    
    P_hat = P_prd;
    
else                         % INS aiding 
    
    % ESKF gain: K[k]
    K = P_prd * Cd' * invQR(Cd * P_prd * Cd' + Rd);
    IKC = eye(15) - K * Cd;
    
    % Estimation error: eps[k]
    eps_pos   = y_pos - p_ins;
    eps_theta = ssa(y_ahrs - theta_ins);  % smallest signed angle   
    
    if (nargin == 10)
        eps = [eps_pos; eps_theta];
    else
        eps_vel = y_vel - v_ins;
        eps = [eps_pos; eps_vel; eps_theta];        
    end
    
    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * eps;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';
    
    % INS reset: x_ins[k]
	p_ins = p_ins + delta_x_hat(1:3);           % Reset INS position
	v_ins = v_ins + delta_x_hat(4:6);           % Reset INS velocity
	b_acc_ins = b_acc_ins + delta_x_hat(7:9);   % Reset INS ACC bias
	theta_ins = theta_ins + delta_x_hat(10:12); % Reset INS attitude
	b_ars_ins = b_ars_ins + delta_x_hat(13:15); % Reset INS ARS bias   
    
end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: x_ins[k+1]
a_ins = R * f_ins + g_n;                             % Linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;           % Exact discretization
v_ins = v_ins + h * a_ins;                           % Exact discretization
theta_ins = theta_ins + h * T * w_ins;               % Euler's method

x_ins = [p_ins; v_ins; b_acc_ins; theta_ins; b_ars_ins];

 
end