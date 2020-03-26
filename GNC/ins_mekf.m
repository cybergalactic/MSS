function [x_ins, P_prd] = ins_mekf(...
   x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos, y_vel)
% Error-state (indirect) feedback EKF for INS aided by position
% measurements y_pos. The velocity aiding signal y_vel is optionally. 
% The EKF is implemented as a feedback filter using reset. 
% Attitude is parametrized using unit quaternions and Gibbs vector is used
% to represent the error states in the MEKF.  
%
% See ExINS_MEKF for how to implement the Kalman filter loop:
%
% if (new slow measurement)
%
%   ins_mekf(...
%   x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos, y_vel)
%
% else (no measurement)
%
%   ins__mekf(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref)
%
% end
%
% Error-state-space model (15 states):
%   delta_x[k+1] = f(delta_x[k], u[k], w[k])
%     delta_y[k] = h(delta_x[k], u[k]) + varepsilon[k]
%
% Inputs: 
%    x_ins[k], INS state vector [p_ins; v_ins, b_acc_ins, q_ins, b_ars_ins] 
%    P_prd[k],  15 x 15 covariance matrix
%    mu, lattitude [rad]
%    h sampling time [s]
%    Qd, Rd, Kalman filter process and measurement covariance matrices
%    f_imu[k], w_imu[k], m_imu[k] fast IMU measurements (spec. force, ARS, mag.)
%    m_ref, magentometer NED reference vector
%           (m_ref = R' m_imu can be computed for phi = theta = psi = 0)
%    y_pos[k], slow position measurement (x, y, z) used for aiding
%    y_vel[k], slow (optionally) velocity measurement used for aiding
%
% Outputs: 
%    x_ins[k+1] state estimate, equal to the INS state vector (16 states)
%    P_prd[k+1] Kalman fiter covariance matrix for error states (15 x 15)
%
% Author:    Thor I. Fossen
% Date:      26 March 2020
% Revisions: 

% Bias time constants (user specified)
T_acc = 1000; 
T_ars = 500; 

%% KF states and matrices
p_ins = x_ins(1:3);          % INS states
v_ins = x_ins(4:6);
b_acc_ins = x_ins(7:9);
q_ins = x_ins(10:13);
b_ars_ins = x_ins(14:16);

% Gravity vector
g = gravity(mu);            % WGS-84 gravity model
g_n = [0 0 g]';    

% Constants 
O3 = zeros(3,3);
I3 = eye(3);

% Attitude matrices
R = Rquat(q_ins);
T = Tquat(q_ins);

% Bias compensated IMU measurements
f_ins = f_imu - b_acc_ins;
w_ins = w_imu - b_ars_ins;

% Magnetometer projections pi_b and pi_n
pi_b = sqrt(m_imu' * m_imu) * eye(3) - m_imu * m_imu';
pi_n = sqrt(m_ref' * m_ref) * eye(3) - m_ref * m_ref';
m_b  = pi_b * m_imu;
m_n  = pi_n * m_ref;

% Reference vectors v10 and v20
v10 = [0 0 1]';                 % gravity vector
v20 = m_n / sqrt( m_n' * m_n);  % magnetic field 
v1 = -f_ins/g;             
v1 = v1 / sqrt( v1' * v1 );
v2 = m_b / sqrt( m_n' * m_n);

% Discrte-time KF matrices
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
if (nargin == 10)             % no aiding
    
    P_hat = P_prd;
    
else                         % INS aiding 
    
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
	p_ins = p_ins + delta_x_hat(1:3);	         % position
	v_ins = v_ins + delta_x_hat(4:6);			 % velocity
	b_acc_ins = b_acc_ins + delta_x_hat(7:9);    % acc bias
	b_ars_ins = b_ars_ins + delta_x_hat(13:15);  % ars bias
     
    q_ins = quatprod(q_ins, delta_q_hat);   % Schur product    
    q_ins = q_ins / sqrt(q_ins' * q_ins);   % normalization       
    
end

% Predictor: P_prd[k+1]
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% INS propagation: x_ins[k+1]
p_ins = p_ins + h * v_ins;
v_ins = v_ins + h * (R * f_ins + g_n);
q_ins = q_ins + h * T * w_ins;
q_ins = q_ins / sqrt(q_ins' * q_ins); % normalization

x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];

end