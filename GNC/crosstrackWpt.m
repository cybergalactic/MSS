function y_e = crosstrackWpt(x2, y2, x1, y1, x, y)
% y_e = crosstrackWpt(x2, y2, x1, y1, x, y, flag) computes the cross-track
% error y_e for a craft at position (x, y) when the path is a straight line
% between two waypoints (x1, y1) and (x2, y2).
%
%
% Input:    (x1,y1) and (x2,y2), waypoints  expressed in NED
%           (x,y), craft North-East positions 
%           flag (optionally) = 1 (plot geometry)
%
% Outputs:  y_e cross-track error 
%
% See also [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag), 
% which computes the coordinate origin (x_p, y_p) of the path-tangential
% reference frame, and the cross-track error y_e.
%  
% Author:    Thor I. Fossen
% Date:      6 Oct. 2020
% Revisions: 

% path-tanegntial angle
pi_p = atan2(y2-y1, x2-x1);

% cross-track error
y_e = -(x-x1) * sin(pi_p) + (y-y1) * cos(pi_p);

