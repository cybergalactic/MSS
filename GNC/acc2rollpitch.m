function [phi, theta] = acc2rollpitch(f)
% acc2rollpitch is compatible with MATLAB and GNU Octave (www.octave.org). 
% The function [phi, theta] = acc2rollpitch(f) computes the static 
% roll-pitch angles phi and theta from 3-axis specific force measurements 
% f = [fx, fy, fz] for multiple sets of measurements (n x 3) or one 
% measurement (1 x 3) or (3 x 1).
%
% Author:    Thor I. Fossen
% Date:      2020-03-20
% Revisions:
%   2024-07-28 : Modified to accept time-series of specific force.

% Input validation and reshaping if necessary
[n, m] = size(f);
if n == 3 && m == 1
    f = f'; % Transpose to a row vector if input is 3x1
    n = 1; % Set n to 1 to indicate a single measurement
elseif n == 1 && m == 3
    % No changes needed for 1x3 input
elseif m ~= 3
    error('Input must have 3 columns representing [fx, fy, fz]');
end

% Initialize output vectors
phi = zeros(n, 1);
theta = zeros(n, 1);

% Compute roll (phi) and pitch (theta) angles for each row
for i = 1:n
    phi(i) = atan2(f(i, 2), f(i, 3));  
    theta(i) = atan2(f(i, 1), sqrt(f(i, 2)^2 + f(i, 3)^2)); 
end

end