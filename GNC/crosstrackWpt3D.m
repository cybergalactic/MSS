function p_e = crosstrackWpt3D(p2, p1, p)
% crosstrackWpt3D is compatible with MATLAB and GNU Octave (www.octave.org). 
% p_e = crosstrackWpt3D(p2, p1, p) computes the 3-D tracking errors
% 
%   p_e = R_y' * R_z' * (p - p1)
%  
% expressed in a path-tangential reference frame P for a craft loacted at 
% the NED position p = [x y z] when the path is a straight line between 
% the waypoints p1 = [x1 y1 z1] and p2 = [x2 y2 z2] expressed in NED.
%
% Input:    p1 = [x1 y1 z1] and p2 = [x2 y2 z2], waypoints  expressed in NED
%           p  = [x y z], craft NED position vector
%
% Outputs:  p_e = [x_e y_e z_e]', along-, cross-, and vertical-track errors
%           expressed in the path-tangential frame P
%
% Author:    Thor I. Fossen
% Date:      2023-10-03
% Revisions: 

% Tracking errors expressed in NED
e_x = p2(1) - p1(1);
e_y = p2(2) - p1(2);
e_z = p2(3) - p1(3);

% Azimuth and elevation angles with respect to the NED frame
pi_h = atan2( e_y, e_x );
pi_v = atan2( -e_z, sqrt( e_x^2 + e_y^2) );

% Rotation matrices
R_z = [ cos(pi_h)   -sin(pi_h)   0
       -sin(pi_h)    cos(pi_h)   0
        0             0          1 ];

R_y = [ cos(pi_v)   0       sin(pi_v) 
        0           1       0
       -sin(pi_v)   0       cos(pi_v) ];

% Tracking errors expressed in path-tangential reference frame P
p_e = R_y' * R_z' * (p - p1);



