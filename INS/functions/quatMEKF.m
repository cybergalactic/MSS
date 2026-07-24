function [quat, b_ars, P_prd] = quatMEKF( ...
    quat, b_ars, P_prd, h, Qd, Rd, T_ars, m_ref, imu_meas)
% quatMEKF is compatible with MATLAB and GNU Octave (www.octave.org).
% This function computes the updated unit quaternion q[k+1], representing 
% the orientation between the BODY and NED frames, as well as the bias 
% b_ars[k+1] of the attitude rate sensor (ARS) in a high-performance 
% Multiplicative Extended Kalman Filter (MEKF) for quaternion attitude 
% estimation (Markley and Crassidis, 2014). The MEKF uses high-rate inertial
% measurements from a 9-DOF inertial measurement unit (IMU). The function 
% supports both predictor mode (IMU only) and corrector mode (IMU with 
% aiding measurements).
%
%   % Predictor (6-DOF IMU only)
%      [quat, b_ars, P_prd] = quatMEKF(quat, b_ars, P_prd, h, Qd, Rd, m_ref, ... 
%          [f_imu', w_imu'])
%   % Corrector (6-DOF IMU + compass aiding)
%      [quat, b_ars, P_prd] = quatMEKF(quat, b_ars, P_prd, h, Qd, Rd, m_ref, ... 
%          [f_imu', w_imu', psi])
%   % Corrector (9-DOF IMU + magnetometer aiding)
%      [quat, b_ars, P_prd] = quatMEKF(quat, b_ars, P_prd, h, Qd, Rd, m_ref, ... 
%          [f_imu', w_imu', m_imu'])
% 
% The quaternion product injection term is implemented using two reference 
% vectors v01 (gravity) and v02 (magnetic field or compass).
% 
% Inputs:   
%   quat[k]   - 4x1 vector representing the current quaternion quat[k] 
%               estimate
%   b_ars[k]  - 3x1 vector representing the current bias of the attitude
%               rate sensor (ARS)
%   P_prd[k]  - 6x6 error-state covariance matrix 
%   h         - Sampling time for the observer update 
%   Qd        - 6x6 process covariance matrix
%   Rd        - Measurement covariance matrix:
%                   6x6 for gravity and magnetometer aiding,
%                   4x4 for gravity and compass aiding.
%   T_ars     - Angular rate bias time constant in seconds.
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
%   Crassidis, J. L., F. L. Markley and Y. Cheng (2007). Survey of Nonlinear 
%       Attitude Estimation Methods. Journal of Guidance, Control and Dynamics 
%       30(1), 12-28.
%
%   Markley, F. L. and J. L. Crassidis (2014). Fundamentals of Spacecraft 
%       Attitude Determination and Control, Volume 33 of Space Technology 
%       Library. Springer-Verlag, New York.
%
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 3rd Edition, Wiley.
%
% Author:    Thor I. Fossen
% Date:      2025-11-05
% Revisions: 

% Transposed unit quaternion rotation matrix: R_transposed[k]
R = Rquat(quat);
R_transposed = R';

% Constants 
O3 = zeros(3,3);
I3 = eye(3);

% High-rate IMU measurements: w_imu[k]
w_imu = imu_meas(4:6)';

% Measurements (low rate)
if numel(imu_meas) == 6  % No magnetic field/compass measurements

    P_hat = P_prd;

else

    if numel(imu_meas) == 9                  % 9‑DOF case: [f_imu' w_imu' m_imu']
        f_imu = imu_meas(1:3)';
        v1  = f_imu / norm(f_imu);
        v01 = [0 0 -1]';

        m_imu = imu_meas(7:9)';
        v2  = m_imu / norm(m_imu);
        v02 = m_ref / norm(m_ref);

        % Measurement matrix
        Cd  = [ Smtrx(R_transposed*v01) O3    % Gravity measurement vector
                Smtrx(R_transposed*v02) O3 ]; % Magnetic field measurement vector

        % Innovation vector
        delta_y = [ v1 - R_transposed*v01
                    v2 - R_transposed*v02];

    elseif numel(imu_meas) == 7               % 7‑DOF case: [f_imu' w_imu' psi]
        f_imu = imu_meas(1:3)';
        v01 = [0 0 -1]';
        v1  = f_imu / norm(f_imu);

        psi = imu_meas(7);
        c_psi = [0, R(3,2), R(3,3)]' / ( R(3,2)^2 + R(3,3)^2 );

        % Measurement matrix
        Cd = [ Smtrx(R_transposed*v01) O3    % Gravity measurement vector
               c_psi'          zeros(1,3) ]; % Compass measurement

        % Innovation vector
        delta_y = [v1 - R_transposed * v01
                   ssa(psi - atan2(R(2,1), R(1,1)))];
    end

    % Kalman gain: K[k])
    K = P_prd * Cd' * invQR(Cd * P_prd * Cd' + Rd);
    IKC = eye(size(P_prd)) - K * Cd;

    % Corrector: delta_x_hat[k] and P_hat[k]
    delta_x_hat = K * delta_y;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';

    % Convert 2 x Gibbs vector to an error quaternion: delta_q_hat[k]
    delta_a = delta_x_hat(1:3);
    delta_q_hat = 1 / sqrt(4 + delta_a' * delta_a) * [2; delta_a];

    % ARS bias reset
    b_ars = b_ars + delta_x_hat(4:6);

    % Multiplicative quaternion update (maintains unit norm)
    quat = quatprod(quat, delta_q_hat);         % Quaternion error injection
    quat = quat / norm(quat);                   % Normalization

end

% ==============================================================================
% Discrete-time ESKF state and process noise matrices Ad and Ed
% ==============================================================================
A = [ -Smtrx(w_imu-b_ars) -I3
       O3                 -1/T_ars * I3 ];

Ad = expm(A * h);

Ed = h *[ -I3 O3
           O3 I3 ];

% ==============================================================================
% Predictor: P_prd[k+1]
% ==============================================================================
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% quat[k+1] is computed using the matrix exponential, which serves as the 
% exponential map for matrix Lie groups, ensuring an exact discretization 
% of the quaternion differential equation: 
%    quat _dot = Tquat(w_imu-b_ars) * quat 
% You can replace the build-in Matlab function expm.m with the custom-made 
% MSS function expm_squaresPade.m for this computation.
quat  = expm(Tquat(w_imu-b_ars) * h) * quat ;   
quat = quat / norm(quat);                        % Normalization

end
