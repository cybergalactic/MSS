function [x_ins, P_prd] = ins_heave(x_ins, P_prd, h, Qd, Rd, f_imu, ...
    phi, theta, p_0, p)
% The function implements an error-state (indirect) feedback Kalman filter 
% (ESKF) in heave. The Inertial Navigation System (INS) is aided by pressure 
% measurements. Usage scenarios are detailed in SIMaidedINSheave, demonstrating 
% the implementation of the Kalman filter loop using the corrector-predictor 
% representation:
%
%   - With new slow pressure measurements:
%       [x_ins,P_prd] = ins_heave(...
%           x_ins, P_prd, h, Qd, Rd, f_imu_z, phi, theta, p_0, p)
%
%   - Without new pressure measurements (no aiding):
%       [x_ins,P_prd] = ins_heave(x_ins, P_prd, h, Qd, Rd, f_imu_z, phi, theta)
%
% This function models the INS errors in heave, including the down position, 
% down velocity, and acceleration bias errors:
%
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k]
%
% Inputs:
%   x_ins[k]   : 3x1 INS state vector at step k, includes down position, down 
%                velocity, and accelerometer biases
%   P_prd[k]   : 3x3 covariance matrix of the prediction step.
%   h          : Sampling time in seconds
%   Qd, Rd     : Process and measurement noise covariance matrices
%   f_imu[k]   : 3 x 1 IMU acceleration in the body frame
%   phi[k]     : Roll angle from AHRS in radians
%   theta[k]   : Pitch angle from AHRS in radians
%   p_0        : Air pressure at the surface in Pa
%   p[k]       : Measured pressure in Pa (aiding measurement)
%
% Outputs:
%   x_prd      : 3x1 Predicted state vector
%   P_prd      : 3x3 Predicted covariance matrix

persistent ins; % Persistent data structure 'ins'

% Initialization of INS parameters. Compute the INS parameters only once
% to avoid that the ESKF repeats the computation in the loop at each time step
if isempty(ins)
    % Constants
    ins.g = 9.81; % Gravity in m/s^2
    ins.rho = 1025; % Density of water in kg/m^3
    ins.T_acc = 1000; % Acceleration bias time constant in seconds

    % Discrete-time ESKF matrices
    ins.A = [ 0 1  0   
              0 0 -1           
              0 0 -1/ins.T_acc ];

    ins.Ad = expm_taylor(ins.A * h);
    ins.Cd = [1 0 0];
    ins.Ed = h * [ 0 0
                  -1 0
                   0 1 ];
end

%% ESKF states
z_ins = x_ins(1);       % Vertical position (NED)
v_z_ins = x_ins(2);     % Vertical velocity (NED)
b_acc_ins = x_ins(3);   % Accelerometer bias in the z direction

% Vertical specific force expressed in NED 
f_D = -f_imu(1) *  cos(phi) * sin(theta) ... 
    + f_imu(2) * sin(phi) ...
    + f_imu(3) * cos(phi) * cos(theta);

% Bias compensated vertical acceleration
a_z_ins = f_D - b_acc_ins + ins.g; 

%% Kalman filter algorithm       
if (nargin == 8)  
    % No aiding    
    P_hat = P_prd;    
else  
    % ESKF gain: K[k]
    K = P_prd * ins.Cd' / (ins.Cd * P_prd * ins.Cd' + Rd);
    IKC = eye(3) - K * ins.Cd;

    % Estimation error: eps_z[k]
    z_meas = (p - p_0) / (ins.rho * ins.g); % p = p_0 + rho * g * z
    eps_z = z_meas - z_ins;

    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * eps_z;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';

    % INS reset: x_ins[k]
	z_ins = z_ins + delta_x_hat(1);           % Reset INS position
	v_z_ins = v_z_ins + delta_x_hat(2);       % Reset INS velocity
	b_acc_ins = b_acc_ins + delta_x_hat(3);   % Reset INS ACC bias	
end

% Predictor: P_prd[k+1]
P_prd = ins.Ad * P_hat * ins.Ad' + ins.Ed * Qd * ins.Ed';

% INS propagation: x_ins[k+1]
z_ins = z_ins + h * v_z_ins + h^2/2 * a_z_ins;  % Exact discretization
v_z_ins = v_z_ins + h * a_z_ins;                % Exact discretization

x_ins = [z_ins v_z_ins b_acc_ins]';

end