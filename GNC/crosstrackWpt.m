function y_e = crosstrackWpt(x2, y2, x1, y1, x, y)
% crosstrack is compatible with MATLAB and GNU Octave (www.octave.org). 
% y_e = crosstrackWpt(x2, y2, x1, y1, x, y, flag) computes the cross-track
% error y_e expressed in NED for a craft located at the North-East position 
% (x, y) when the path is a straight line from waypoint (x1,y1) to 
% waypoint (x2,y2).
%
% Input:    (x1,y1) and (x2,y2): Waypoints  expressed in NED
%           (x,y): Craft North-East positions 
%
% Outputs:  y_e: Cross-track error expressed in NED
%
% See also [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag), 
% which computes the coordinate origin (x_p, y_p) of the path-tangential
% reference frame, and the cross-track error y_e.
%  
% Author:    Thor I. Fossen
% Date:      2020-10-06
% Revisions: 
%   None

% Path-tanegntial angle with respect to the North axis
pi_h = atan2(y2-y1, x2-x1);

% Cross-track error expressed in NED
y_e = -(x-x1) * sin(pi_h) + (y-y1) * cos(pi_h);

