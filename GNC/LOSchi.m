function [chi_d, omega_chi_d,y_e] = LOSchi(x,y,Delta,R_switch,wpt,U,K_f,h)
% [chi_d, omega_chi_d,y_e] = LOSchi(x,y,Delta,R_switch,wpt,U,K_f,h)
% LOSpsi computes thedesired course angle when the path is straight lines 
% going through the waypoints  (wpt.pos.x, wpt.pos.y). The desired course
% angle chi_d and course rate d/dt chi_d = omega_chi_d used by course 
% autopilot systems are computed using the roportional LOS guidance law:
%
%  chi_d = pi_h - atan( Kp * y_e ),    Kp = 1/Delta  
%
%  omega_chi_d = -Kp * U * sin( chi - pi_h ) / ( 1 + (Kp * y_e)^2 )
%
% where pi_h is the path-tangential (azimuth) angle with respect to the North 
% axis and y_e is the cross-track error. The observer/filter for the desired 
% yaw angle psi_d is
%
%    d/dt psi_f = r_d + K_f * ssa( psi_d - psi_f ) )
%
% where the desired course rate omega_d = d(psi_d)/dt is computed by
%
%    d/dt y_e = -U * y_e / sqrt( Delta^2 + y_e^2 )
%    omega_chi_d = -Kp * dy_e/dt / ( 1 + (Kp * y_e)^2 )  
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
%   U: speed, vehicle cruise speed or time-varying measurement (m/s)
%   K_f: observer gain for desired yaw angle (typically 0.1-0-5)
%   h: sampling time (s)
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%      dist = sqrt(  (wpt.pos.x(k+1)-wpt.pos.x(k))^2 
%                  + (wpt.pos.y(k+1)- wpt.pos.y(k))^2 );
%
% Outputs:  
%    chi_d:       desired course angle (rad)
%    omega_chi_d: desired course rate (rad/s)
%
% For heading control use the functions LOSpsi.m and ILOSpsi.m.
%
% Ref. T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
% Motion Control. 2nd. Edition, Wiley
%
% Author:    Thor I. Fossen
% Date:      2 June 2021
% Revisions: 18 June 2021 - added output omega_chi_d
%            17 Oct 2022  - added filter/observer for chi_d to avoid steps 

persistent k;      % active waypoint index (initialized by: clear LOSchi)
persistent xk yk;  % active waypoint (xk, yk) corresponding to integer k
persistent chi_f;

%% Initialization of (xk, yk) and (xk_next, yk_next)
if isempty(k)   
  
    % check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end
    
    % check input parameters
    if (R_switch < 0); error("R_switch must be larger than zero"); end
    if (Delta < 0); error("Delta must be larger than zero"); end    
    
    k = 1;              % set first waypoint as the active waypoint
    xk = wpt.pos.x(k);
    yk = wpt.pos.y(k);   
    chi_f = 0;          % filtered chi_d
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

%% Print active waypoint 
fprintf('Active waypoint:\n')
fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);

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
end

% LOS guidance law
chi_d = chi_f;
chi_ref = pi_h - atan( y_e/Delta );

% desired course rate: omega_chi_d[k]
Dy_e = -U * y_e / sqrt( Delta^2 + y_e^2 );        % d/dt y_e
omega_chi_d = -(1/Delta) * Dy_e / ( 1 + (y_e/Delta)^2 );   

% observer: chi_f[k+1]
chi_f = chi_f + h * (omega_chi_d + K_f * ssa( chi_ref - chi_f ) );

end

