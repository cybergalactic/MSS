function [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, imu_meas)
% quatObserver is compatible with MATLAB and GNU Octave (www.octave.org).
% This function computes the updated unit quaternion q[k+1], representing 
% the orientation between the BODY and NED frames, as well as the bias 
% b_ars[k+1] of the attitude rate sensor (ARS) in a high-performance 
% nonlinear observer (Grip et al 2013). The observer uses high-rate inertial
% measurements from a 9-DOF inertial measurement unit (IMU). The function 
% can be called either as a corrector (with new measurements) or as a 
% predictor (without new measurements). Additionally, the magnetometer 
% can operate at a slower rate (typically 100 Hz) compared to the high-rate 
% specific force and ARS measurements (typically 1000 Hz). The function 
% should be called from a synchronized 'Main Loop' running at twice the 
% high-rate frequency (e.g., sampling time h = 1/2000 corresponding to 
% 2000 Hz) to satisfy the Nyquist sampling theorem.
%
%   Corrector with new f_imu, w_imu and m_imu measurements:
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu', m_imu'])
%   Corrector with new f_imu and w_imu measurements:
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref, ... 
%          [f_imu', w_imu'])
%
%   Predictor (no measurements):
%      [quat, b_ars] = quatObserver(quat, b_ars, h, Ki, k1, k2, m_ref)
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
%   h         - Time step (sampling period) for the observer update
%   Ki        - 3x3 diagonal integral gain matrix for bias correction
%   k1        - Gain for the feedback term associated with the specific 
%               force measurement
%   k2        - Gain for the feedback term associated with the magnetic 
%               field measurement
%   m_ref     - 3x1 vector representing the reference magnetic field vector
%               in the NED frame. The reference signal m_ref = R^n_b * m_imu 
%               can be computed during initial calibration, see 
%               staticRollPitchYaw.m
% imu_meas[k] = [f_imu', w_imu', m_imu'] is 1x9 vector with components
%    [fx, fy, fz, wx, wy, wz, mx, my, mz]. More specific, f_imu[k] is a
%    3x1 vector representing the IMU specific force measurements, w_imu[k] 
%    is a 3x1 vector representing the IMU angular velocity measurements, and 
%    m_imu[k] is a 3x1 vector representing the IMU magnetic field measurements.
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
if nargin == 8

    % Specific force and angular velocity IMU measurements
    f_imu = imu_meas(1:3);
    w_imu = imu_meas(4:6);

    v01 = [0 0 1]'; % Normalized gravity reference vector (NED)
    v1 = -f_imu / norm(f_imu); % Normalized specific force (forward-starboard-down)

    % Nonlinear injection terms
    sigma1 = k1 * cross(v1, R_transposed * v01);  % Specific force
    sigma2 = 0; % Magnetic field
    
    % If new magnetometer measurements, [f_imu w_imu m_imu] has 9 columns
    [~, N] = size(imu_meas);
    if N == 9  

        % Magnetic field IMU measurements
        m_imu = imu_meas(7:9);

        v02 = m_ref / norm(m_ref); % Normalized magnetic field reference vector (NED)
        v2 =  m_imu / norm(m_imu); % Normalized magnetic field (forward-starboard-down)
        
        % Make the observer robust to situations where the magnetic field vector is
        % nearly aligned with gravity using a threhold value for vector projection
        if abs(dot(v2, v1)) <= 0.9   
            % Projection of the vector x2 onto the plane orthogonal to x1 is
            % x2_projected = x2 - dot(x2, x1) * x1 whenever norm(x) = 1
            v2 = v2 - dot(v2, v1) * v1;
            v02 = v02 - dot(v02, v01) * v01;

            % Normalize the projected vectors
            v2 = v2 / norm(v2);
            v02 = v02 / norm(v02);
        end

        % Nonlinear injection term
        sigma2 = k2 * cross(v2, R_transposed * v02); % Magnetic field
    end

% Nonlinear injection term
sigma = sigma1 + sigma2; % Specific force + magnetic field

end

% State propagation: quat[k+1] and b_ars[k+1]
b_ars = b_ars - h * Ki * sigma; % Attitude rate sensor (ARS) bias
quat = expm( Tquat(w_imu - b_ars + sigma) * h ) * quat; % Unit quaternion

% Normalization
quat = quat / norm(quat); 

end
