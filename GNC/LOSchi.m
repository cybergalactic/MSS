function [chi_d, omega_chi_d] = LOSchi(x,y,Delta,R_switch,wpt,U,chi)
% [chi_d, omega_chi_d] = LOSchi(x,y,Delta,R_switch,wpt,U,chi) computes the
% desired course angle when the path is straight lines going through the 
% waypoints  (wpt.pos.x, wpt.pos.y). The desired course angle chi_d and 
% course rate d/dt \chi_d = omega_chi_d (optionally) used by course 
% autopilot systems are computed using the roportional LOS guidance law:
%
%  chi_d = pi_p - atan( Kp * y_e ),    Kp = 1/Delta  
%
%  omega_chi_d = -Kp * U * sin( chi - pi_p ) / ( (Kp * y_e)^2 + 1);
%
% where pi_p is the path-tangential angle with respect to the North axis
% and y_e is the cross-track error expressed in NED. The function can be
% called according to:
%
%  chi_d = LOSchi(x,y,Delta,R_switch,wpt)
%  [chi_d, omega_chi_d] = LOSchi(x,y,Delta,R_switch,wpt,U,chi)
%
% Initialization:
%   The active waypoint (xk, yk) where k = 1,2,...,n is a persistent
%   integer should be initialized to the first waypoint, k = 1, using 
%   >> clear LOSchi
%
% Inputs:   
%   (x,y): craft North-East positions (m)
%   Delta: positive look-ahead distance (m)
%   R_switch: go to next waypoint when the along-track distance x_e 
%             is less than R_switch (m)
%   wpt.pos.x = [x1, x2,..,xn]': array of waypoints expressed in NED (m)
%   wpt.pos.y = [y1, y2,..,yn]': array of waypoints expressed in NED (m)
%   U: speed (m/s), only needed for computation of omega_chi_d
%   chi: course angle (rad), only needed for computation of omega_chi_d
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%      dist = sqrt(  (wpt.pos.x(k+1)-wpt.pos.x(k))^2 
%                  + (wpt.pos.y(k+1)- wpt.pos.y(k))^2 );
%
% Outputs:  
%    chi_d:       desired course angle (rad)
%    omegs_chi_d: desired course rate (rad/s)
%
% For integral LOS (course control) use ILOSchi.m. 
% For heading control use the functions LOSpsi.m and ILOSpsi.m.
%
% Ref. T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
% Motion Control. 2nd. Edition, Wiley
%
% Author:    Thor I. Fossen
% Date:      2 June 2021
% Revisions: 18 June 2021 - added optional output omega_chi_d

persistent k;   % active waypoint index (initialized by: clear LOSchi)
persistent xk;  % active waypoint (xk, yk) corresponding to integer k
persistent yk;

%% Initialization of (xk, yk) and (xk_next, yk_next)
if isempty(k)   
  
    % check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end
    
    % check input parameters
    if (R_switch < 0)
        error("R_switch must be larger than zero");
    end
    if (Delta < 0)
        error("Delta must be larger than zero");
    end    
    
    k = 1;              % set first waypoint as the active waypoint
    xk = wpt.pos.x(k);
    yk = wpt.pos.y(k);     
end

%% Read next waypoint (xk_next, yk_next) from wpt.pos 
n = length(wpt.pos.x);
if k < n                        % if there are more waypoints, read next one 
    xk_next = wpt.pos.x(k+1);  
    yk_next = wpt.pos.y(k+1);    
else                            % else, use the last one in the array
    xk_next = wpt.pos.x(end);
    yk_next = wpt.pos.y(end); 
end

%% Print active waypoint 
fprintf('Active waypoint:\n')
fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);

%% Compute the desired course angle
pi_p = atan2(yk_next-yk, xk_next-xk);  % path-tangential angle w.r.t. to North

% along-track and cross-track errors (x_e, y_e) expressed in NED
x_e =  (x-xk) * cos(pi_p) + (y-yk) * sin(pi_p);
y_e = -(x-xk) * sin(pi_p) + (y-yk) * cos(pi_p);

% if the next waypoint satisfy the switching criterion, k = k + 1
d = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
if ( (d - x_e < R_switch) && (k < n) )
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
end

% LOS guidance law
Kp = 1/Delta;
chi_d = pi_p - atan( Kp * y_e );

% kinematic differential equations
Dy_e = U * sin( chi - pi_p );

% Course rate (optionally)
if (nargin == 7)
    omega_chi_d = -Kp * Dy_e / ( (Kp * y_e)^2 + 1);
end

end

