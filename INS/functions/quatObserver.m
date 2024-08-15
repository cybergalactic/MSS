function [quat, b_ars] = quatObserver(...
    quat, b_ars, h, Ki, k1, k2, g, f_imu, w_imu, m_imu, m_ref)
% quatObserver is compatible with MATLAB and GNU Octave (www.octave.org).
% This function computes the updated unit quaternion q[k+1], representing the 
% orientation between the BODY and NED frames, and the bias b_ars[k+1] of the 
% attitude rate sensor (ARS) in a nonlinear observer using inertial 
% measurements from a 9-DOF inertial measurement unit (IMU). The function 
% can be called as either a corrector (with new measurements) or a predictor 
% (without new measurements):
%
%   Corrector (new IMU measurements)
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, g, ...
%         f_imu, w_imu, m_imu, m_ref)
%
%   Predictor (no measurements)
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, g)
% 
% The injection term is implemented using two reference vectors
%   sigma = k1 * v1 x R'(q) * v01 + k2 * v2 x R'(q) * v02
%
% Continuous-time observer (Grip et al. 2013)(Fossen 2021, Eqs. 14.48-14.50)
%   quat_dot = T(q) * (w_imu - b_ars - sigma)
%   b_ars_dot = -Ki * sigma
%
% Discrete-time observer using exact discretization (Fossen 2021), where the
% injection term, sigma = 0 (predictor) and sigma ~= 0 (corrector):
%   quat = expm(Tquat(w_imu - b_ars + sigma) * h) * quat
%   quat = quat / sqrt(quat' * quat)
% 
% Inputs:   
%   quat[k]   - 4x1 vector representing the current quaternion quat[k] estimate
%   b_ars[k]  - 3x1 vector representing the current bias of the attitude
%               rate sensor (ARS)
%   f_imu[k]  - 3x1 vector representing the IMU specific force measurement 
%               f_imu[k] in the BODY frame 
%   w_imu[k]  - 3x1 vector representing the IMU angular velocity measurement
%               in the BODY frame 
%   m_imu[k]  - 3x1 vector representing the IMU magnetic field measurement 
%               in the BODY frame 
%   m_ref     - 3x1 vector representing the reference magnetic field vector
%               in the NED frame. The reference signal m_ref = R^n_b * m_imu 
%               can be computed during initial calibration by aligning the 
%               IMU such that phi = theta = psi = 0. Then m_ref = m_imu.
%   Ki        - 3x3 diagonal integral gain matrix for bias correction
%   k1        - Gain for the feedback term associated with the specific 
%               force measurement
%   k2        - Gain for the feedback term associated with the magnetic 
%               field measurement
%   g         - Acceleration due to gravity, typically 9.81 m/s^2
%   h         - Time step (sampling period) for the observer update
%
% Outputs:  
%   quat[k+1]  - 4x1 vector representing the updated quaternion estimate
%   b_ars[k+1] - 3x1 vector representing the updated ARS bias estimate
%
% The function implements a quaternion-based nonlinear observer for attitude
% estimation by Grip et al. (2013), which guarantees USGES stability. 
% The observer updates the quaternion based on specific force and magnetic 
% field measurements while also correcting for ARS bias. The observer employs 
% a feedback mechanism that uses the cross products between the measured vectors 
% and their reference counterparts to compute an injection term, which is used 
% to update the bias and quaternion. USGES stability guarantees robustness 
% to bounded disturbances. 
%
% References:
%   H. F. Grip, T. I. Fossen, T. A. Johansen, and A. Saberi (2013). Nonlinear 
%       Observer for GNSS-Aided Inertial Navigation with Quaternion-Based 
%       Attitude Estimation. American Control Conference (ACC), Washington DC, 
%       USA, IEEE Xplore, pp. 272-279. doi.org/10.1109/ACC.2013.6579849
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 2nd Edition, Wiley.
%
% Author:    Thor I. Fossen
% Date:      2024-08-15
% Revisions: 

% Transposed unit quaternion rotation matrix
R_transposed = Rquat(quat);

% Injection term: sigma[k]
sigma = 0; % Observer is a predictor by default

% Observer acts as a corrector when there are new measurements
if nargin == 11 
    
    % NED normalized reference vectors 
    v01 = [0 0 1]'; % Gravity reference vector
    v02 = m_ref / sqrt( m_ref' * m_ref ); % Magnetic field reference vector 

    % BODY-fixed normalized IMU measurements
    v1 = -f_imu / g; % Specific force measurement 
    v1 = v1 / sqrt( v1' * v1 ); 
    v2  = m_imu / sqrt( m_imu' * m_imu ); % Magnetic field measurement

    % Nonlinear injection term
    sigma = k1 * cross(v1, R_transposed * v01) ...
          + k2 * cross(v2, R_transposed * v02);
end

% State propagation: quat[k+1] and b_ars[k+1]
b_ars = b_ars - h * Ki * sigma; % Attitude rate sensor (ARS) bias
quat = expm( Tquat(w_imu - b_ars + sigma) * h ) * quat; % Unit quaternion

% Normalization
quat = quat / sqrt(quat' * quat); 

end
