function [phi, theta, psi] = staticRollPitchYaw(f_imu, m_imu)
% staticRollPitchYaw is compatible with MATLAB and GNU Octave (www.octave.org). 
% This function computes the static roll-pitch-yaw angles (phi, theta, psi)
% from 3-axis specific force and magnetometer measurements expressed in the 
% BODY frame. If only specific force is used as input the function returns
% the static roll and pitch angles. 
% 
% A useful application is computation of the magentic field reference vector 
%
%    m_ref = Rzyx(phi,theta,psi) * m_imu
%
% expressed in NED, which is used by quatObserver.m to compute the 
% attitude of a moving rigid-body. This assumes that the IMU is at rest 
% initially when phi, theta, and psi are computed. 
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
%   psi   : Nx1 vector of yaw angles in radians.
%
% Author:    Thor I. Fossen
% Date:      2024-08-17
% Revisions:

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

% Initialize output vectors
phi = zeros(rows1, 1);
theta = zeros(rows1, 1);
psi = zeros(rows1, 1);

% Compute roll (phi), pitch (theta), and yaw (psi) angles for each row
for i = 1:rows1
    % Static roll and pitch angles, see Fossen (2021, Eqs. (14.34-(14.35)).
    % Calculate roll angle using atan(fy/fz) since atan2(fy, fz) fails when 
    % both arguments can have both signs. This approach avoids 180 deg roll 
    % angles when the physical angle is 0 deg.
    phi(i) = atan(f(i, 2) / f(i, 3));  

    % Correct the angle based on the quadrant
    if f(i, 3) > 0  % fz > 0, meaning potential Quadrant 2 or 3
        if f(i, 2) > 0
            phi(i) = phi(i) - pi;  % Quadrant 2: 90 to 180 degrees
        else
            phi(i) = phi(i) + pi;  % Quadrant 3: -90 to -180 degrees
        end
    end

    % Calculate pitch angle using atan2. This works for both signs of fx 
    % since sqrt(fy^2 + fz^2) > 0 (always positive).
    theta(i) = atan2(f(i, 1), sqrt(f(i, 2)^2 + f(i, 3)^2)); 

    % Compute yaw angle (psi) only if magnetometer data is provided
    if nargin == 2
        % Tilt-compensated magnetometer readings, see Fossen (2021, Eq.(14.27)).
        % [mx, my, mz]' = R_y(theta) * R_x(phi) * [m_imu_x, m_imu_y, m_imu_z]'
        mx = m_imu(i, 1) * cos(theta(i)) ...
            + m_imu(i, 2) * sin(phi(i)) * sin(theta(i)) ...
            + m_imu(i, 3) * cos(phi(i)) * sin(theta(i));
        my = m_imu(i, 2) * cos(phi(i)) ...
            - m_imu(i, 3) * sin(phi(i));
        
        % Static yaw angle
        psi(i) = atan2(my, mx);
    end
end

end