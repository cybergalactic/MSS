function [quat, b_ars] = quatObserver( ...
    quat, b_ars, h, Ki, k1, k2, m_ref, imu_meas, coningSculling)
% quatObserver is compatible with MATLAB and GNU Octave (www.octave.org).
% This function computes the updated unit quaternion q[k+1], representing 
% the orientation between the BODY and NED frames, as well as the bias 
% b_ars[k+1] of the attitude rate sensor (ARS) in a high-performance 
% nonlinear observer (Grip et al 2013). The observer uses high-rate inertial
% measurements from a 9-DOF inertial measurement unit (IMU). The function 
% can be called either as a corrector (with new measurements) or as a 
% predictor (without new measurements). 
%
%   9-DOF measurements:
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu', m_imu'])
%   7-DOF measurements:
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu', psi])
%   6-DOF measurements (no magnetometer/compass measurements):
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu'])
% 
% The injection term is implemented using two reference vectors
%   sigma = k1 * v1 x R'(quat) * v01 + k2 * v2 x R'(quat) * v02
%
% Continuous-time observer (Grip et al. 2013)(Fossen 2021, Eqs. 14.48-14.50)
%   quat_dot = Tquat(w_imu - b_ars + sigma) * quat
%   b_ars_dot = -Ki * sigma
%
% Discrete-time observer using the matrix exponential, which serves as the 
% exponential map for matrix Lie groups, ensuring an exact discretization 
% of the quaternion differential equation: 
%   quat = expm( Tquat(w_imu - b_ars + sigma) * h ) * quat
%   quat = quat / sqrt(quat' * quat)
% 
% Inputs:   
%   quat[k]   - 4x1 vector representing the current quaternion quat[k] 
%               estimate
%   b_ars[k]  - 3x1 vector representing the current bias of the attitude
%               rate sensor (ARS)
%   h         - Sampling time for the observer update 
%   Ki        - 3x3 diagonal integral gain matrix for bias estimation
%   k1        - Gain for the injections term associated with the specific 
%               force measurement vector
%   k2        - Gain for the injection term associated with the magnetic 
%               field measurement vector
%   m_ref     - 3x1 vector representing the reference magnetic field vector
%               expressed in NED. The reference signal m_ref = R^n_b * m_imu 
%               can be computed during initial calibration, see 
%               staticRollPitchYaw.m
%  imu_meas[k] = [f_imu', w_imu', m_imu'] is 1x9 vector with components
%               [fx, fy, fz, wx, wy, wz, mx, my, mz]. More specific, 
%               f_imu[k] is a 3x1 vector representing the IMU specific force 
%               measurements, w_imu[k] is a 3x1 vector representing the IMU 
%               angular velocity measurements, and m_imu[k] is a 3x1 vector 
%               representing the IMU magnetic field measurements. The IMU 
%               axes are assumed to be oriented forward-starboard-down.
%               For compass measurements [f_imu', w_imu', psi] is of dimension 7. 
% coningSculling (OPTIONAL) 1, to apply coning and sculling compensation.
%                0, is the default value (no compensation). Coning is the 
%                error from integrating small rotational movements, while 
%                sculling is the error from combining rotational movements 
%                with linear accelerations during integration. The effects
%                are small when using high-rate measurements (1000 Hz).
%                Compensation is important for slower update rates.
%
% Outputs:  
%   quat[k+1]  - 4x1 vector representing the updated unit quaternion estimate
%   b_ars[k+1] - 3x1 vector representing the updated ARS bias estimate
%
% The function implements the quaternion-based nonlinear observer for 
% attitude estimation by Grip et al. (2013). The observer updates the 
% quaternion based on specific force and magnetic field measurements while 
% also correcting for ARS bias. The observer employs a feedback mechanism 
% that uses the cross products between the measured vectors and their 
% reference counterparts to compute an injection term, which is used to 
% update the bias and quaternion. USGES stability guarantees robustness to
% bounded disturbances. 
%
% References:
%   H. F. Grip, T. I. Fossen, T. A. Johansen, and A. Saberi (2013). 
%       Nonlinear Observer for GNSS-Aided Inertial Navigation with 
%       Quaternion-Based Attitude Estimation. American Control Conference, 
%       Washington DC, USA, IEEE Xplore, pp. 272-279. 
%       doi.org/10.1109/ACC.2013.6579849
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 2nd Edition, Wiley.
%
% Author:    Thor I. Fossen
% Date:      2024-08-20
% Revisions: 
%   2025-06-17 Modified to accept compass measurements.

if nargin == 8
    coningSculling = 0; % Default, no compensation of coning and sculling
end

% Transposed unit quaternion rotation matrix: R_transposed[k]
R_transposed = Rquat(quat)';

% High-rate IMU specific force and ARS measurements: f_imu[k] and w_imu[k]
f_imu = imu_meas(1:3)';
w_imu = imu_meas(4:6)';

% Nonlinear injection terms (low rate)
if length(imu_meas) == 9
    % IMU measuremenst with 3-axis magnetometer measurement: [f_imu' w_imu' m_imu'] 
    m_imu = imu_meas(7:9)'; % Magnetic field IMU measurements: m_imu[k] in BODY
    v01 = [0 0 -1]'; % Normalized gravity reference vector in NED (measuring -g at rest)
    v1 = f_imu / norm(f_imu); % Normalized specific force measurement
    v02 = m_ref / norm(m_ref); % Normalized magnetic field reference vector in NED
    v2 =  m_imu / norm(m_imu); % Normalized magnetic field measurement in BODY
    sigma1 = k1 * cross(v1, R_transposed * v01);
    sigma2 = k2 * cross(v2, R_transposed * v02);
    sigma = sigma1 + sigma2;

elseif length(imu_meas) == 7 
    % IMU measurements with scalar compass meausrement: [f_imu' w_imu' psi]  
    psi = imu_meas(7);  % Compass measurement: psi[k]
    v01 = [0 0 -1]'; % Normalized gravity reference vector in NED (measuring -g at rest)
    v1 = f_imu / norm(f_imu); % Normalized specific force measurement
    v02 = [1; 0; 0];  % Normalized heading reference vector in NED
    v2 = [cos(psi); -sin(psi); 0];  % Normalized heading measurement vector in BODY
    sigma1 = k1 * cross(v1, R_transposed * v01);
    sigma2 = k2 * cross(v2, R_transposed * v02);
    sigma = sigma1 + sigma2;

else
    % No correction
    sigma = zeros(3,1);
end

% State propagation: quat[k+1] is computed using the matrix exponential, 
% which serves as the exponential map for matrix Lie groups, ensuring an 
% exact discretization of the quaternion differential equation: 
%   quat_dot = Tquat(w_imu - b_ars + sigma) * quat
% You can replace the build-in Matlab function expm.m with the custom-made 
% MSS functions expm_taylor.m or expm_squaresPade.m for this computation.

% Angular velocity with bias compensation 'b_ars' and injection term 'sigma'
w_estimated = w_imu - b_ars + sigma; 

if coningSculling == 0
    % No coning and sculling compensation 
    quat = expm( Tquat(w_estimated) * h ) * quat; % Quaternion propagation
    quat = quat / norm(quat);  % Normalization
else
    % Midtpoint method for compensating coning and sculling effects 
    quat_midpoint = expm( Tquat(w_estimated) * h/2) * quat;
    quat_midpoint = quat_midpoint / norm(quat_midpoint);
    quat = expm( Tquat(w_estimated) * h/2 ) * quat_midpoint;
    quat = quat / norm(quat);  % Normalization
end

% State propagation: b_ars[k+1]
b_ars = b_ars - h * Ki * sigma; % Attitude rate sensor (ARS) bias

end
