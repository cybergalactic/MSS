function [x_ins, P_prd] = ins_mekf(...
   x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos, y_vel)
% ins_mekf is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter 
% (ESKF) specifically for Inertial Navigation Systems (INS) that are 
% aided by magnetometer and positional data. Attitude is parametrized using 
% the 4-parameter unit quaternion representation and the Gibbs vector in 
% the Multiplicative Error State Kalman Filter (MEKF) formulation, thus 
% avoiding gimbal lock. 
%
% Usage scenarios are detailed in examples SIMaidedINSquat and ExINS_MEKF, 
% demonstrating the implementation of the Kalman filter loop using the 
% corrector-predictor representation:
%
%   - With new slow measurements:
%      [x_ins,P_prd] = ins_mekf(...
%         x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos)
%      [x_ins,P_prd] = ins_mekf(...
%         x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos, y_vel)
%
%   - Without new measurements (no aiding):
%      [x_ins,P_prd] = ins_mekf(...
%         x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref)
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
%   m_imu[k] : Magnetic field measurements from the IMU.
%   m_ref    : Magentometer NED reference vector
%             (m_ref = R' m_imu can be computed for phi = theta = psi = 0)
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
% Date: 2020-04-26
% Revisions: 
%   2020-11-29: Bugfix - removed magnetometer projection algorithms.
%   2021-12-13: Improved numerical stability by replacing Euler's method 
%               with exact discretization in the INS state propagation.

% Bias time constants (user specified)
T_acc = 1000; 
T_ars = 500; 

%% ESKF states and matrices
p_ins = x_ins(1:3);          % INS states
v_ins = x_ins(4:6);
b_acc_ins = x_ins(7:9);
q_ins = x_ins(10:13);
b_ars_ins = x_ins(14:16);

% Gravity vector
g = gravity(mu);             % WGS-84 gravity model
g_n = [0 0 g]';    

% Constants 
O3 = zeros(3,3);
I3 = eye(3);

% Unit quaternion rotation matrix
R = Rquat(q_ins);

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Normalized gravity vectors
v10 = [0 0 1]';             % NED
v1 = -f_ins/g;              % BODY
v1 = v1 / sqrt( v1' * v1 );

% Normalized magnetic field vectors
v20 = m_ref / sqrt( m_ref' * m_ref );    % NED
v2  = m_imu / sqrt( m_ref' * m_ref );    % BODY
 
% Discrte-time ESKF matrices
A = [ O3 I3  O3            O3              O3
      O3 O3 -R            -R*Smtrx(f_ins)  O3
      O3 O3 -(1/T_acc)*I3  O3              O3
      O3 O3  O3           -Smtrx(w_ins)   -I3
      O3 O3  O3            O3             -(1/T_ars)*I3 ];
   
Ad = eye(15) + h * A + 0.5 * (h * A)^2;

if (nargin == 11)
    Cd = [ I3 O3 O3 O3 O3               % NED positions 
           O3 O3 O3 Smtrx(R'*v10) O3    % Gravity           
           O3 O3 O3 Smtrx(R'*v20) O3 ]; % Magnetic field
else
    Cd = [ I3 O3 O3 O3 O3               % NED positions 
           O3 I3 O3 O3 O3               % NED velocities 
           O3 O3 O3 Smtrx(R'*v10) O3    % Gravity          
           O3 O3 O3 Smtrx(R'*v20) O3 ]; % Magnetic field
end

Ed = h *[  O3 O3    O3 O3
          -R  O3    O3 O3
           O3 I3    O3 O3
           O3 O3   -I3 O3
           O3 O3    O3 I3  ];

%% Kalman filter algorithm       
if (nargin == 10)             % No aiding
    
    P_hat = P_prd;
    
else                          % INS aiding 
    
    % KF gain: K[k]
    K = P_prd * Cd' * inv(Cd * P_prd * Cd' + Rd); 
    IKC = eye(15) - K * Cd;
    
    % Estimation error: eps[k]
    eps_pos = y_pos - p_ins;
    eps_g   = v1 - R' * v10; 
    eps_mag = v2 - R' * v20;    
    
    if (nargin == 11)
        eps = [eps_pos; eps_g; eps_mag];
    else
        eps_vel = y_vel - v_ins;
        eps = [eps_pos; eps_vel; eps_g; eps_mag];        
    end
    
    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * eps;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';
    
	% Error quaternion (2 x Gibbs vector): delta_q_hat[k]
	delta_a = delta_x_hat(10:12);
	delta_q_hat = 1/sqrt(4 + delta_a' * delta_a) * [2 delta_a']';
    
    % INS reset: x_ins[k]
	p_ins = p_ins + delta_x_hat(1:3);	         % Reset position
	v_ins = v_ins + delta_x_hat(4:6);			 % Reset velocity
	b_acc_ins = b_acc_ins + delta_x_hat(7:9);    % Reset ACC bias
	b_ars_ins = b_ars_ins + delta_x_hat(13:15);  % Reset ARS bias
     
    q_ins = quatprod(q_ins, delta_q_hat);   % Schur product    
    q_ins = q_ins / sqrt(q_ins' * q_ins);   % Normalization       
    
end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: x_ins[k+1]
a_ins = R * f_ins + g_n;                     % Linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;   % Exact discretization
v_ins = v_ins + h * a_ins;                   % Exact discretization
q_ins = expm( Tquat(w_ins) * h ) * q_ins;    % Exact discretization
% q_ins = q_ins + h * Tquat(q_ins) * w_ins;  % Euler's method (alternative)
q_ins = q_ins / sqrt(q_ins' * q_ins);        % Normalization

x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];

end