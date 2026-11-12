function [quat, b_ars] = quatObserver(quat,b_ars,h,Ki,k1,k2,m_ref,imu_meas)
% This function computes the updated unit quaternion q[k+1], representing 
% the orientation between the BODY and NED frames, as well as the bias 
% b_ars[k+1] of the attitude rate sensor (ARS) in a high-performance 
% nonlinear observer (Mahony, Hamel and Pflimlin, 2008) (Grip et al 2013). 
% The function  supports both predictor mode (IMU only) and corrector mode 
% (IMU with aiding measurements).
%
%   % Predictor (6-DOF IMU only)
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu'])
%   % Corrector (6-DOF IMU + compass aiding)
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu', psi])
%   % Corrector (9-DOF IMU with magnetometer aiding)
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu', m_imu'])
% 
% With magnetometer aiding, the injection uses two reference vectors:
%   sigma = k1 * cross(v1, R' * v01) + k2 * cross(v2, R' * v02)
% With scalar compass aiding, the second term is the wrapped yaw error
% about the NED down axis, expressed in BODY.
%
% Continuous-time observer 
%   quat_dot = Tquat(w_imu - b_ars + sigma) * quat
%   b_ars_dot = -Ki * sigma
%
% Discrete-time observer with the corrected angular rate held constant
% over each sample interval:
%   quat = expm( Tquat(w_imu - b_ars + sigma) * h ) * quat
%   quat = quat / sqrt(quat' * quat)
%   b_ars = b_ars - h * Ki * sigma
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
%   k2        - Gain for the magnetic field vector or scalar compass heading
%   m_ref     - 3x1 vector representing the reference magnetic field vector
%               expressed in NED. The reference signal m_ref = R^n_b * m_imu 
%               can be computed during initial calibration, see 
%               staticRollPitchYaw.m
%  imu_meas[k] - Measurement vector:
%                  [f_imu', w_imu', m_imu'] (9×1) accel, gyro, mag
%                  [f_imu', w_imu', psi]    (7×1) accel, gyro, compass yaw
%                  [f_imu', w_imu']         (6×1) accel, gyro only
%               f_imu[k] : 3x1 High-rate IMU specific force measurements.
%               w_imu[k] : 3x1 High-rate IMU angular velocity measurements. 
%               m_imu[k] : 3x1 Low-rate IMU magnetic field measurements. 
%               psi[k]   : 1x1 Low-rate compass measurement.
% 
%               The IMU axes are assumed to be oriented forward-starboard-down.
%
% Outputs:  
%   quat[k+1]  - 4x1 vector representing the updated unit quaternion estimate
%   b_ars[k+1] - 3x1 vector representing the updated ARS bias estimate
%
% References:
%   H. F. Grip, T. I. Fossen, T. A. Johansen, and A. Saberi (2013). 
%       Nonlinear Observer for GNSS-Aided Inertial Navigation with 
%       Quaternion-Based Attitude Estimation. American Control Conference, 
%       Washington DC, USA, IEEE Xplore, pp. 272-279. 
%       doi.org/10.1109/ACC.2013.6579849
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 3rd Edition, Wiley.
%
% Author:    Thor I. Fossen
% Date:      2024-08-20
% Revisions: 
%   2025-06-17 Modified to accept compass measurements.
%   2026-11-06 Inject the compass heading error about the NED down axis.

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
    % IMU measurements with scalar compass measurement: [f_imu' w_imu' psi]  
    psi = imu_meas(7);  % Compass measurement: psi[k]
    v01 = [0 0 -1]'; % Normalized gravity reference vector in NED (measuring -g at rest)
    v1 = f_imu / norm(f_imu); % Normalized specific force measurement
    
    % Inject the wrapped yaw error about the NED down axis expressed in BODY.
    % The injection vanishes at the correct attitude wherever yaw is defined.
    % The former BODY-forward comparison could introduce a spurious correction
    % at nonzero pitch even with exact heading measurements.
    R = R_transposed';                       % BODY to NED
    psi_hat = atan2(R(2,1), R(1,1));         % Estimated yaw angle
    sigma1 = k1 * cross(v1, R_transposed * v01);
    sigma2 = k2 * ssa(psi - psi_hat) * (R_transposed * [0; 0; 1]);
    sigma = sigma1 + sigma2;

else
    % No correction
    sigma = zeros(3,1);
end

% State propagation uses the matrix exponential with the corrected angular
% rate held constant over the sample interval. MSS alternatives to expm are
% expm_taylor.m and expm_squaresPade.m.

% Angular velocity with bias compensation 'b_ars' and injection term 'sigma'
w_estimated = w_imu - b_ars + sigma; 

% Propagate the quaternion using the exponential map.
quat = expm( Tquat(w_estimated) * h ) * quat; 
quat = quat / norm(quat);  % Normalization

% State propagation: b_ars[k+1]
b_ars = b_ars - h * Ki * sigma; % Attitude rate sensor (ARS) bias

end
