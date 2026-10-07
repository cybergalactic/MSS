function [phi, theta, psi] = staticRollPitchYaw(f_imu, m_imu)
% This function computes static roll, pitch, and magnetic heading
% (phi, theta, psi) from 3-axis specific force and magnetometer measurements
% expressed in the BODY frame. If only specific force is used as input, the
% function returns the static roll and pitch angles.
% The BODY axes are forward-starboard-down, and the navigation axes are NED.
%
% For positive-east magnetic declination delta, the true NED yaw angle and
% magnetic field reference vector are
%
%    psi_true = ssa(psi + delta)
%    m_ref = Rzyx(phi,theta,psi_true) * m_imu
%
% The NED reference vector can be used by quatObserver.m to compute the
% attitude of a moving rigid body. This assumes that the IMU is initially
% at rest when phi, theta, and psi are computed.
%
% Syntax:
%   [phi, theta, psi] = staticRollPitchYaw(f_imu, m_imu)
%   [phi, theta] = staticRollPitchYaw(f_imu)
%
% Inputs:
%   f_imu : A matrix of size Nx3, where each row contains the specific force 
%           measurements [fx, fy, fz] expressed in BODY.
%   m_imu : (Optional) A matrix of size Nx3, where each row contains the 
%           magnetometer measurements [mx, my, mz] expressed in BODY.
%
% Outputs:
%   phi   : NX1 vector of roll angles in radians.
%   theta : Nx1 vector of pitch angles in radians.
%   psi   : Nx1 vector of magnetic yaw angles in radians.
%
% Author:    Thor I. Fossen
% Date:      2024-08-17
% Revisions:
%   2026-10-07 : Corrected IMU references and FSD/NED heading sign (E. Krizman)

% Input validation
narginchk(1, 2); % Ensure at least one input and at most two inputs are provided

[rows1, cols1] = size(f_imu);
if cols1 ~= 3
    error('f_imu should have 3 columns corresponding to [fx, fy, fz].');
end

if nargin == 2
    [rows2, cols2] = size(m_imu);
    if cols2 ~= 3
        error('m_imu should have 3 columns corresponding to [mx, my, mz].');
    end
    if rows1 ~= rows2
        error('f_imu and m_imu must have the same number of rows.');
    end
end

% Compute static roll and pitch angles
[phi, theta] = acc2rollpitch(f_imu);

% Initialize magnetic heading vector
psi = zeros(rows1, 1);

% Compute magnetic heading only if magnetometer data is provided
if nargin == 2
    for i = 1:rows1
        % Tilt-compensated magnetometer readings, see Fossen (2027,
        % Eqs. (14.15)-(14.16)).
        % [hx, hy, hz]' = R_y(theta) * R_x(phi) * [m_imu_x, m_imu_y, m_imu_z]'
        hx = m_imu(i, 1) * cos(theta(i)) ...
            + m_imu(i, 2) * sin(phi(i)) * sin(theta(i)) ...
            + m_imu(i, 3) * cos(phi(i)) * sin(theta(i));
        hy = m_imu(i, 2) * cos(phi(i)) ...
            - m_imu(i, 3) * sin(phi(i));
        
        % For FSD BODY axes, the starboard component has the opposite sign
        % of positive clockwise heading relative to magnetic north.
        psi(i) = atan2(-hy, hx);
    end
end

end
