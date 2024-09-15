function [phi, theta] = acc2rollpitch(f_imu, b_acc)
% acc2rollpitch is compatible with MATLAB and GNU Octave (www.octave.org). 
% The function [phi, theta] = acc2rollpitch(f_imu, b_acc) computes the static 
% roll-pitch angles phi and theta from 3-axis specific force measurements 
% f_imu = [fx, fy, fz] for multiple sets of measurements (n x 3) or one 
% measurement (1 x 3) or (3 x 1). The 3x1 vector b_acc is an optional
% acceleration bias compensation term satisfying f = f_imu - b_acc. For the NED
% reference frame, the IMU should measure f = [0 0 -g]' at rest when levelled.
%
% Author:    Thor I. Fossen
% Date:      2020-03-20
% Revisions:

% Input validation and reshaping if necessary
[n, m] = size(f_imu);
if n == 3 && m == 1
    f_imu = f_imu'; % Transpose to a row vector if f_imu is 3x1
    n = 1; % Set n to 1 to indicate a single measurement
elseif m ~= 3
    error('Input must have 3 columns representing [fx, fy, fz]');
end

if nargin == 1
    b_acc = zeros(1,3);
else
    b_acc = b_acc(:)'; % Ensure that b_acc is a row vector
end

% Initialize output vectors
phi = zeros(n, 1);
theta = zeros(n, 1);

% Compute roll (phi) and pitch (theta) angles for each row
for i = 1:n
    
    % Bias-compensated specific force
    f = f_imu(i, :) - b_acc;

    % Calculate roll and pitch angles 
    phi(i) = atan(f(2) / f(3));  
    theta(i) = atan(f(1) / sqrt(f(2)^2 + f(3)^2));  

end

end