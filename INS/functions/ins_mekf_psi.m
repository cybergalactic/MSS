function [x_ins, P_prd] = ins_mekf_psi(...
    x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, f_imu, w_imu, y_psi, y_pos, y_vel)
% ins_mekf_psi is compatible with MATLAB and GNU Octave (www.octave.org).
% The function implements an error-state (indirect) feedback Kalman filter 
% (ESKF) specifically for Inertial Navigation Systems (INS) that are 
% aided by compass and positional data. Attitude is parameterized using the 
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
%          x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, f_imu, w_imu, y_psi, y_pos)
%       [x_ins,P_prd] = ins_mekf_psi(...
%          x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, f_imu, w_imu, y_psi, y_pos, y_vel)
%
%   - % Without new position measurements:
%       [x_ins,P_prd] = ins_mekf_psi(...
%          x_ins, P_prd, mu, h, Qd, Rd, T_acc, T_ars, f_imu, w_imu, y_psi)
%
% This function models the INS errors in a 15-dimensional state space, 
% including position, velocity, biases, and attitude errors:
%
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k]
%
% Inputs:
%   x_ins[k] : INS state vector at step k, includes position, velocity, 
%              accelerometer biases, attitude (quaternion), and gyro biases.
%   P_prd[k] : 15x15 covariance matrix of the prediction step.
%   mu       : Latitude in radians, used to calculate Earth's gravity vector.
%   h        : Sampling time in seconds.
%   Qd, Rd   : Process and measurement noise covariance matrices for the 
%              Kalman filter.
%   T_acc    : Acceleration bias time constant in seconds.
%   T_ars    : Angular rate bias time constant in seconds.
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
%   T. I. Fossen (2027). "Handbook of Marine Craft Hydrodynamics and Motion 
%   Control," 3rd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2020-04-26
% Revisions: 
%   2021-12-21: Improved numerical accuracy by replacing Euler's method 
%               with exact discretization in the INS state propagation.
%   2022-08-30: Use atan2 instead of atan to avoid jumps in the formula:
%               eps_psi = ssa( y_psi - atan2(u_y, u_x) ); 
%   2024-09-09 : Redesign for slower compass measurements
%   2025-11-19 : Add new aiding methods/measurements options

% ==============================================================================
% ESKF signals
% ==============================================================================
p_ins = x_ins(1:3);          % INS states
v_ins = x_ins(4:6);
b_acc_ins = x_ins(7:9);
q_ins = x_ins(10:13);
b_ars_ins = x_ins(14:16);
if Rd.pseudoFlag, zn_int_ins = x_ins(17); end % Additional sea level state

% WGS-84 gravity vector expressed in NED
g_n = [0 0 gravity(mu)]';   

% Matrix constants
O3 = zeros(3,3);
I3 = eye(3);

% Unit quaternion rotation matrix
Rq = Rquat(q_ins);

% ==============================================================================
% Discrete-time ESKF measurement matrix Cd
% ==============================================================================
Cd = [];
delta_y = [];
Rd.mtrx = [];
       
% 1) Position aiding (if y_pos exists)
if nargin >= 12
    Cd = [Cd; I3 O3 O3 O3 O3];
    delta_y = [delta_y; y_pos - p_ins];
    Rd.mtrx = blkdiag(Rd.mtrx, Rd.position);
end

% 2) Velocity aiding (if y_vel exists)
if nargin == 13
    Cd = [Cd; O3 I3 O3 O3 O3];
    delta_y = [delta_y; y_vel - v_ins];
    Rd.mtrx = blkdiag(Rd.mtrx, Rd.velocity);
end

% 3) Linearization of the compass measurement equation if y_psi is received
if ~isempty(y_psi)
    c_psi = [0, Rq(3,2), Rq(3,3)] / ( Rq(3,2)^2 + Rq(3,3)^2 );
    Cd = [Cd; zeros(1,9) c_psi zeros(1,3)];
    delta_y = [delta_y; ssa(y_psi - atan2(Rq(2,1), Rq(1,1)))];
    Rd.mtrx = blkdiag(Rd.mtrx, Rd.compass);
end

% 4) Gravity-reference aiding when translational acceleration a_nM is negligible
if Rd.applicationFlag
    v01 = [0 0 -1]'; % NED reference vector (measuring -g at rest)
    f_gravity = f_imu - b_acc_ins; % Bias-compensated force 
    v1_gravity  = f_gravity / norm(f_gravity); % Specific force measurement
    Cd = [Cd; O3 O3 O3 Smtrx(Rq'*v01) O3];
    delta_y = [delta_y; v1_gravity - Rq'*v01];
    Rd.mtrx = blkdiag(Rd.mtrx, Rd.gravityRefVector);
end

% 5) Sea level pseudo-measurement (integral of zn is 0)
if Rd.pseudoFlag
    Cd = [Cd zeros(size(Cd,1),1); zeros(1,15) 1];
    delta_y = [delta_y; -zn_int_ins];
    Rd.mtrx = blkdiag(Rd.mtrx, Rd.pseudoMeasSeaLevel);
end

% ==============================================================================
%% ESKF      
% ==============================================================================
if isempty(delta_y)           % No aiding measurement
    P_hat = P_prd;  
else                          % Aiding 
    % KF gain: K[k]
    K = P_prd * Cd' / (Cd * P_prd * Cd' + Rd.mtrx); 
    IKC = eye(size(P_prd)) - K * Cd;
    
    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * delta_y;
    P_hat = IKC * P_prd * IKC' + K * Rd.mtrx * K';
    
    % INS reset: x_ins[k]
	p_ins = p_ins + delta_x_hat(1:3);	          % Reset position
	v_ins = v_ins + delta_x_hat(4:6);			  % Reset velocity
	b_acc_ins = b_acc_ins + delta_x_hat(7:9);     % Reset ACC bias
	b_ars_ins = b_ars_ins + delta_x_hat(13:15);   % Reset ARS bias

    % Convert 2 x Gibbs vector to an error quaternion
    delta_a = delta_x_hat(10:12);
    delta_q_hat = [2; delta_a] / sqrt(4 + delta_a' * delta_a);

    if Rd.pseudoFlag
 	   zn_int_ins = zn_int_ins + delta_x_hat(16); % Reset integral of vertical position
    end

    q_ins = quatprod(q_ins, delta_q_hat);        % Quaternion error injection    
    q_ins = q_ins / norm(q_ins);                 % Normalization           
end

% ==============================================================================
% Bias-corrected IMU measurements
% ==============================================================================
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% ==============================================================================
% Discrete-time ESKF state and process noise matrices Ad and Ed
% ==============================================================================
Rq = Rquat(q_ins); % Recompute Rq using updated q_ins
 
A = [ ...
    O3 I3  O3            O3               O3
    O3 O3 -Rq           -Rq*Smtrx(f_ins)  O3
    O3 O3 -(1/T_acc)*I3  O3               O3
    O3 O3  O3           -Smtrx(w_ins)    -I3
    O3 O3  O3            O3              -(1/T_ars)*I3 ];

if Rd.pseudoFlag
    A = [ A zeros(15,1)
        zeros(1,16) ];
end

Ad = expm(A * h); 

Ed = h *[ ...
    O3 O3    O3 O3
   -Rq O3    O3 O3
    O3 I3    O3 O3
    O3 O3   -I3 O3
    O3 O3    O3 I3  ];
if Rd.pseudoFlag
    Ed = blkdiag(Ed, 1);
end

% ==============================================================================
% Predictor: P_prd[k+1]
% ==============================================================================
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% ==============================================================================
% INS propagation: p_ins[k] and v_ins[k]
% ==============================================================================
a_ins = Rq * f_ins + g_n;                    % Linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;   % Exact discretization
v_ins = v_ins + h * a_ins;                   % Exact discretization

% q_ins[k+1] is computed using the matrix exponential, which serves as the 
% exponential map for matrix Lie groups, ensuring an exact discretization 
% of the quaternion differential equation: 
%    q_ins_dot = Tquat(w_ins) * q_ins
% You can replace the build-in Matlab function expm.m with the custom-made 
% MSS function expm_squaresPade.m for this computation.
q_ins = expm(Tquat(w_ins) * h) * q_ins;      % Exponential map
q_ins = q_ins / norm(q_ins);                 % Normalization

% INS state vector: x_ins[k+1]
x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];
if Rd.pseudoFlag
    x_ins = [x_ins; zn_int_ins]; % Additional sea level state
end

end