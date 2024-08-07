function dense_wpt = addIntermediateWaypoints(wpt, multiplier)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org).
% The function addIntermediateWaypoints adds intermediate waypoints along 
% the line segments between given waypoints for better resolution.
%
% Inputs:
%   wpt       - A struct containing fields `pos.x` and `pos.y` that are vectors
%               of x and y coordinates of waypoints.
%   multiplier- Number of segments each original segment should be divided into.
%
% Outputs:
%   dense_wpt - A struct containing fields `pos.x` and `pos.y` with the new,
%               denser set of waypoints.
%
% Author:    Thor I. Fossen
% Date:      2024-04-21
% Revisions:
%   None

% Ensure wpt.pos.x and wpt.pos.y are column vectors
wpt.pos.x = wpt.pos.x(:);
wpt.pos.y = wpt.pos.y(:);
dense_wpt.pos.x = [];
dense_wpt.pos.y = [];

for i = 1:length(wpt.pos.x) - 1
    % Add the current waypoint
    dense_wpt.pos.x(end + 1) = wpt.pos.x(i);
    dense_wpt.pos.y(end + 1) = wpt.pos.y(i);

    % Calculate the number of intermediate points to add
    num_intermediate_points = multiplier - 1;

    % Linear interpolation for intermediate points
    for j = 1:num_intermediate_points
        t = j / multiplier;
        new_x = wpt.pos.x(i) * (1 - t) + wpt.pos.x(i + 1) * t;
        new_y = wpt.pos.y(i) * (1 - t) + wpt.pos.y(i + 1) * t;
        dense_wpt.pos.x(end + 1) = new_x;
        dense_wpt.pos.y(end + 1) = new_y;
    end
end

% Add the last waypoint
dense_wpt.pos.x(end + 1) = wpt.pos.x(end);
dense_wpt.pos.y(end + 1) = wpt.pos.y(end);

end
