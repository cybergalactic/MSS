function [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag)
% crosstrack is compatible with MATLAB and GNU Octave (www.octave.org). 
% [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag) computes 
% the coordinate origin (x_p, y_p) of the path-tangential reference frame.
% The straight-line path goes from reference point (x_ref, y_ref) to the 
% target position (x_p, y_p). The cross-track error y_e for a craft at
% located at (x, y) is the length of the orthogonal vector between the 
% origin (x_p, y_p) and (x, y). 
%
% Inputs:    
%   (x_t, y_t)          : North-East target positions
%   (x_ref, y_ref)      : North-East reference positions
%   (x, y)              : Craft North-East positions
%   flag (optionally)   : 1 to plot the geometry, else no plot
%
% Outputs:  
%   (x_p, y_p)          : Origin of path-tangential frame expressed in NED
%   y_e                 : Cross-track error
%
% Usage:
%   [xp, yp, ye] = crosstrack(150, 300, 0, 0, 80, 30, 1)   % Plot geometry
%   [xp, yp, ye] = crosstrack(150, 300, 0, 0, 80, 30)      % No plot
%
% For GNC applications it is not necessary to compute (xp, yp) implying
% that visualization and matrix inversion can be avoided. Hence, it is
% reccomended to use the computational more efficient functions:
%
%   2-D waypoints: crosstrackWpt(x2, y2, x1, y1, x, y)
%   3-D waypoints: crosstrackWpt3([x2, y2, z2], [x1, y1, z1], [x, y, z]) 
%  
% Author:    Thor I. Fossen
% Date:      2020-07-10
% Revisions: 
%   2020-10-06 - Fixed minor typos
%   2025-05-06 - Replaced matrix inversion method with orthogonal projection

% Path vector (NE)
p1 = [x_ref; y_ref];      % Start waypoint (N, E)
p2 = [x_t; y_t];          % End waypoint (N, E)
p  = [x; y];              % Vehicle position (N, E)
d  = p2 - p1;             % Direction vector (N, E)

% Orthogonal projection of p onto line (p1, p2)
t = dot(p - p1, d) / dot(d, d);
p_p = p1 + t * d;

% Output reference point
x_p = p_p(1);  % North
y_p = p_p(2);  % East

% Path heading (angle from North axis)
pi_h = atan2(d(2), d(1));  % atan2(East, North)

% Cross-track error in NE frame
y_e = -(x - x_p) * sin(pi_h) + (y - y_p) * cos(pi_h);

if (nargin == 7 && flag == 1)
    axis normal;
    figure(gcf)
    plot(y,x,'ro',y_ref,x_ref,'ko',y_t,x_t,'yo',y_p,x_p,'ks',...
        'linewidth',2,'MarkerFaceColor','c','MarkerSize',15)
    hold on
    quiver(y_ref,x_ref,y_t-y_ref,x_t-x_ref,0,'-b','linewidth',2)
    quiver(y_p,x_p,y-y_p,x-x_p,0,'-r','linewidth',2)
    hold off
    grid;
    xlabel('East')
    ylabel('North')
    axis equal;
    legend('Vehicle position', 'Start position','Target position','Origin',...
        'Path','Cross-track error','FontSize',12)
end

end

