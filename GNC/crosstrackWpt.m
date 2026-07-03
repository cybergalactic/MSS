function y_e = crosstrackWpt(x2, y2, x1, y1, x, y)
% crosstrackWpt is compatible with MATLAB and GNU Octave (www.octave.org).
% y_e = crosstrackWpt(x2, y2, x1, y1, x, y) computes the signed cross-track
% error y_e for a craft located at the North-East position (x, y) relative
% to the straight-line path from waypoint (x1, y1) to waypoint (x2, y2).
%
% Inputs:
%   (x1, y1) : Start waypoint expressed in NED
%   (x2, y2) : End waypoint expressed in NED
%   (x, y)   : Craft position expressed in NED
%
% Output:
%   y_e      : Signed cross-track error
%
% See also crosstrack, which additionally computes the orthogonal % projection
% (x_p, y_p) of the craft position onto the path. This point defines the origin 
% of the path-tangential reference frame.
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

% Path-tanegntial angle with respect to the North axis
pi_h = atan2(y2-y1, x2-x1);

% Cross-track error expressed in NED
y_e = -(x-x1) * sin(pi_h) + (y-y1) * cos(pi_h);

