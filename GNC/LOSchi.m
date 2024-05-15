function [chi_ref, y_e] = LOSchi(x, y, Delta_h, R_switch, wpt)
% [chi_ref, y_e]  = LOSchi(x, y, Delta_h, R_switch, wpt)
% LOSchi computes the desired course angle, chi_ref, and cross-track
% error, y_e, for paths that consist of straight lines connecting waypoints 
% (wpt.pos.x, wpt.pos.y). The desired course angle is computed using the 
% adaptive LOS (ALOS) guidance law 
%
%  chi_ref[k] = pi_h - atan( y_e[k] / Delta_h ) 
%
% where pi_h is the path-tangential (azimuth) angle with respect to the 
% North axis and y_e is the cross-track error. To handle steps in the course 
% angle command, chi_ref, due to waypoint switching, the following observer 
% should be used:
%
%  [chi_d, omega_chi_d] = LOSobserver(chi_d omega_chi_d,chi_ref,h,K_f,T_f)
%
% where chi_d is the desired course angle and omega_chi_d isthe desired 
% course rate.  
%
% Initialization:
%   The active waypoint (xk, yk) where k = 1,2,...,n is a persistent
%   integer should be initialized to the first waypoint, k = 1, using 
%   >> clear LOSchi
%
% Inputs:   
%   (x,y):    craft North-East positions (m)
%   Delta_h:  positive look-ahead distance (m)
%   R_switch: go to next waypoint when the along-track distance x_e 
%             is less than R_switch (m)
%   wpt.pos.x = [x1, x2,..,xn]': array of waypoints expressed in NED (m)
%   wpt.pos.y = [y1, y2,..,yn]': array of waypoints expressed in NED (m)
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%     dist = sqrt(   ( wpt.pos.x(k+1) - wpt.pos.x(k) )^2 
%                  + ( wpt.pos.y(k+1) - wpt.pos.y(k) )^2 )
% Outputs:  
%    chi_ref: LOS course angle (rad)
%    y_e:     cross-track error (m)
%
% For heading control use the functions LOSpsi.m and ILOSpsi.m.
% See also: crosstracHermiteLOS.m
%
% Ref. T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
% Motion Control. 2nd. Edition, Wiley
%
% Author:    Thor I. Fossen
% Date:      2021-06-21
% Revisions: 2024-04-01 - Added compability to LOSobserver.m

persistent k;      % active waypoint index (initialized by: clear LOSchi)
persistent xk yk;  % active waypoint (xk, yk) corresponding to integer k

%% Initialization of (xk, yk) and (xk_next, yk_next)
if isempty(k)   
  
    % check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end
    
    % check input parameters
    if (R_switch < 0); error("R_switch must be larger than zero"); end
    if (Delta_h < 0); error("Delta must be larger than zero"); end    
    
    k = 1;              % set first waypoint as the active waypoint
    xk = wpt.pos.x(k);
    yk = wpt.pos.y(k);  
    fprintf('Active waypoints:\n')
    fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);

end

%% Read next waypoint (xk_next, yk_next) from wpt.pos 
n = length(wpt.pos.x);
if k < n                        % if there are more waypoints, read next one 
    xk_next = wpt.pos.x(k+1);  
    yk_next = wpt.pos.y(k+1);    
else                            % else, continue with last bearing
    bearing = atan2((wpt.pos.y(n)-wpt.pos.y(n-1)), (wpt.pos.x(n)-wpt.pos.x(n-1)));
    R = 1e10;
    xk_next = wpt.pos.x(n) + R * cos(bearing);
    yk_next = wpt.pos.y(n) + R * sin(bearing); 
end

%% Compute the desired course angle w.r.t. North
pi_h = atan2( (yk_next-yk), (xk_next-xk) );  

% along-track and cross-track errors (x_e, y_e)
x_e =  (x-xk) * cos(pi_h) + (y-yk) * sin(pi_h);
y_e = -(x-xk) * sin(pi_h) + (y-yk) * cos(pi_h);

% if the next waypoint satisfy the switching criterion, k = k + 1
d = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
if ( (d - x_e < R_switch) && (k < n) )
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
    fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);
end

% LOS guidance law
chi_ref = pi_h - atan( y_e/Delta_h );

end

