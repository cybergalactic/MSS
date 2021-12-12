function [x_ins, P_prd] = ins_euler(...
   x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, psi, y_pos, y_vel)
% Error-state (indirect) feedback KF for INS aided by position and compass
% measurements y_pos and psi. The velocity aiding signal y_vel is optionally. 
% The KF is implemented as a feedback filter using reset. 
% Attitude is parametrized using Euler angles.  
%
% See ExINS_Euler for how to implement the Kalman filter loop:
%
% if (new slow measurement)
%
%   ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, psi, y_pos, y_vel)
%
% else (no measurement)
%
%   ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, psi)
%
% end
%
% Error-state-space model (15 states):
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k]
%
% Inputs: 
%    x_ins[k], INS state vector [p_ins, v_ins, b_acc_ins, q_ins, b_ars_ins] 
%    P_prd[k],  15 x 15 covariance matrix
%    mu, lattitude [rad]
%    h sampling time [s]
%    Qd, Rd, Kalman filter process and measurement covariance matrices
%    f_imu[k], w_imu[k], m_imu[k] fast IMU measurements (spec. force, ARS, mag.)
%    psi[k],   compass measurement
%    y_pos[k], slow position measurement (x, y, z) used for aiding
%    y_vel[k], slow (optionally) velocity measurement used for aiding
%
% Outputs: 
%    x_ins[k+1] state estimate, equal to the INS state vector (15 states)
%    P_prd[k+1] Kalman fiter covariance matrix for error states (15 x 15)
%
% Author:    Thor I. Fossen
% Date:      14 Jan 2021
% Revisions: 12 Dec 2021 - replaced Euler's method with exact discretization

% Bias time constants (user specified)
T_acc = 1000; 
T_ars = 500; 

%% KF states and matrices
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

% Attitude matrices
R = Rzyx(theta_ins(1),theta_ins(2),theta_ins(3));
T = Tzyx(theta_ins(1),theta_ins(2));

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Normalized gravity vectors
v10 = [0 0 1]';             % NED
v1 = -f_ins/g;              % BODY
v1 = v1 / sqrt( v1' * v1 );
 
% Discrete-time KF matrices
A = [ O3 I3  O3            O3              O3
      O3 O3 -R            -R*Smtrx(f_ins)  O3
      O3 O3 -(1/T_acc)*I3  O3              O3
      O3 O3  O3           -Smtrx(w_ins)   -I3
      O3 O3  O3            O3             -(1/T_ars)*I3 ];
   
Ad = eye(15) + h * A + 0.5 * (h * A)^2;

if (nargin == 10)
    Cd = [ I3 O3 O3 O3 O3                % NED positions 
           O3 O3 O3 Smtrx(R'*v10) O3     % Gravity  
           zeros(1,11) 1 zeros(1,3) ];   % Compass
else
    Cd = [ I3 O3 O3 O3 O3                % NED positions 
           O3 I3 O3 O3 O3                % NED velocities 
           O3 O3 O3 Smtrx(R'*v10) O3     % Gravity    
           zeros(1,11) 1 zeros(1,3) ];   % Compass       
end

Ed = h *[  O3 O3    O3 O3
          -R  O3    O3 O3
           O3 I3    O3 O3
           O3 O3   -I3 O3
           O3 O3    O3 I3  ];

%% Kalman filter algorithm       
if (nargin == 9)             % no aiding
    
    P_hat = P_prd;
    
else                         % INS aiding 
    
    % KF gain: K[k]
    K = P_prd * Cd' * inv(Cd * P_prd * Cd' + Rd); 
    IKC = eye(15) - K * Cd;
    
    % Estimation error: eps[k]
    eps_pos = y_pos - p_ins;
    eps_g   = v1 - R' * v10; 
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
	p_ins = p_ins + delta_x_hat(1:3);	         % position
	v_ins = v_ins + delta_x_hat(4:6);			 % velocity
	b_acc_ins = b_acc_ins + delta_x_hat(7:9);    % acc bias
    theta_ins = theta_ins + delta_x_hat(10:12);  % theta
	b_ars_ins = b_ars_ins + delta_x_hat(13:15);  % ars bias
          
end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: x_ins[k+1]
a_ins = R * f_ins + g_n;                        % linear acceleration
p_ins = p_ins + h * v_ins + h^2/2 * a_ins;      % exact discretization
v_ins = v_ins + h * a_ins;                      % exact discretization
theta_ins = theta_ins + h * T * w_ins;          % Euler's method

x_ins = [p_ins; v_ins; b_acc_ins; theta_ins; b_ars_ins];

end