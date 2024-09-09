function [x_ins, P_prd] = ins_mekf_psi(...
    x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos, y_vel)
% ins_mekf_psi is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter 
% (ESKF) specifically for Inertial Navigation Systems (INS) that are 
% aided by compass and positional data. Attitude is parametrized using the 
% 4-parameter unit quaternion representation and the Gibbs vector in the 
% Multiplicative Error State Kalman Filter (MEKF) formulation, thus avoiding
% gimbal lock. 
%
% Usage scenarios are detailed in SIMaidedINSquat.m demonstrating the 
% implementation of the Kalman filter loop using the corrector-predictor 
% representation:
%
%   - With new slow position measurements:
%       [x_ins,P_prd] = ins_mekf_psi(...
%          x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos)
%       [x_ins,P_prd] = ins_mekf_psi(...
%          x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos, y_vel)
%
%   - Without new position measurements (no aiding):
%       [x_ins,P_prd] = ins_mekf_psi(...
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
%   f_imu[k] : High-rate IMU specific force measurements.
%   w_imu[k] : High-rate IMU angular rate measurements. 
%   y_psi[k] : Compass measurement (yaw angle).
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
%   2021-12-21: Improved numerical accuracy by replacing Euler's method 
%               with exact discretization in the INS state propagation.
%   2022-08-30: Use atan2 instead of atan to avoid jumps in the formula:
%               eps_psi = ssa( y_psi - atan2(u_y, u_x) ); 
%   2024-08-31: Sign correction for v10 and using invQR.m instead of inv.m
%   2024-09-09 : Redesign for slower compass measurements

% Bias time constants (user specified)
T_acc = 500; 
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

% NED normalized reference vector
v01 = [0 0 -1]';  % Gravity reference vector, f_z = [0 0 -g]' at rest

% Unit quaternion rotation matrix
R = Rquat(q_ins);

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Discrete-time ESKF matrices
A = [ O3 I3  O3            O3              O3
      O3 O3 -R            -R*Smtrx(f_ins)  O3
      O3 O3 -(1/T_acc)*I3  O3              O3
      O3 O3  O3           -Smtrx(w_ins)   -I3
      O3 O3  O3            O3             -(1/T_ars)*I3 ];
   
% Ad = eye(15) + h * A + 0.5 * (h * A)^2 + ...
Ad = expm_taylor(A * h); 
       
if (nargin == 10)
    Cd = [ I3 O3 O3 O3 O3                 % NED positions (x, y, z)
           O3 O3 O3 Smtrx(R'*v01) O3 ];   % Gravity           
else
    Cd = [ I3 O3 O3 O3 O3                 % NED positions (x, y, z)
           O3 I3 O3 O3 O3                 % NED velocities 
           O3 O3 O3 Smtrx(R'*v01) O3 ];   % Gravity          
end

% Linearization of the compass measurement equation if y_psi is received
if ~isempty(y_psi)
    a = (2/q_ins(1)) * [q_ins(2) q_ins(3) q_ins(4)]'; % 2 x Gibbs vector
    u_y = 2 * ( a(1)*a(2) + 2*a(3) );
    u_x = ( 4 + a(1)^2 - a(2)^2 - a(3)^2 );
    u = u_y / u_x;
    du = 1 / (1 + u^2);
    c_psi = du * (1 / ( 4 + a(1)^2 - a(2)^2 - a(3)^2 )^2 ) * ...
        [ -2*((a(1)^2 + a(3)^2 - 4)*a(2) + a(2)^3 + 4*a(1)*a(3))
        2*((a(2)^2 - a(3)^2 + 4)*a(1) + a(1)^3 + 4*a(2)*a(3))
        4*((a(3)^2 + a(1)*a(2)*a(3) + a(1)^2 - a(2)^2 + 4)) ];

    Cd = [ Cd
        zeros(1,9) c_psi' zeros(1,3)]; % Augment the compass measurement to Cd
end

Ed = h *[  O3 O3    O3 O3
          -R  O3    O3 O3
           O3 I3    O3 O3
           O3 O3   -I3 O3
           O3 O3    O3 I3  ];

%% Kalman filter algorithm       
if (nargin == 9)               % No aiding
    P_hat = P_prd;
else                           % Aiding    
    % KF gain: K[k]
    K = P_prd * Cd' * invQR(Cd * P_prd * Cd' + Rd); 
    IKC = eye(15) - K * Cd;
    
    % BODY-fixed normalized IMU measurements (forward-starboard-down)
    v1 = f_ins / norm(f_ins);  % Specific force measurement
    
    % Estimation errors (injection terms)
    eps_pos = y_pos - p_ins;
    eps_g   = v1 - R' * v01; 
    
    eps_psi = ssa( y_psi - atan2(u_y, u_x) );   
    
    if (nargin == 10)
        eps = [eps_pos; eps_g; eps_psi];
    else
        eps_vel = y_vel - v_ins;
        eps = [eps_pos; eps_vel; eps_g; eps_psi];        
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
     
    q_ins = quatprod(q_ins, delta_q_hat);        % Schur product    
    q_ins = q_ins / norm(q_ins);                 % Normalization          
end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: p_ins[k] and v_ins[k]
a_ins = R * f_ins + g_n;                     % Linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;   % Exact discretization
v_ins = v_ins + h * a_ins;                   % Exact discretization

% q_ins[k+1] is computed using the matrix exponential, which serves as the 
% exponential map for matrix Lie groups, ensuring an exact discretization 
% of the quaternion differential equation: 
%    q_ins_dot = Tquat(w_ins) * q_ins
% You can replace the build-in Matlab function expm.m with the custom-made 
% MSS function expm_squaresPade.m for this computation.
q_ins = expm( Tquat(w_ins) * h ) * q_ins;    % Exact discretization
q_ins = q_ins / norm(q_ins);                 % Normalization

x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];

end