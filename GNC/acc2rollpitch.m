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
%   2024-08-17 : Made robust for four-quadrant angle solution

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

end

end