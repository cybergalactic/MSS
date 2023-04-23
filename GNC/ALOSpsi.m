function [psi_d, r_d, y_e] = ALOSpsi(x,y,Delta,gamma,h,R_switch,wpt,U,K_f)
% [psi_d, r_d, y_e] = ALOSpsi(x,y,Delta,gamma,h,R_switch,wpt,U) 
% ALOSpsi computes the desired heading angle psi_d, desired yaw rate r_d 
% and cross-track error y_e when the path is  straight lines going through 
% the waypoints (wpt.pos.x, wpt.pos.y). The desired heading angle computed 
% using the adaptive LOS guidance law by Fossen (2023).
%
%    psi_d = pi_h - beta_hat - atan( Kp * y_e ),  Kp = 1/Delta
%
%    d/dt beta_hat = gamma * Delta * y_e / sqrt( Delta^2 + y_e^2 )
%
% where pi_h is the path-tangential (azimuth) angle with respect to the 
% North axis and y_e is the cross-track error. The observer/filter for the 
% desired yaw angle psi_d is
%
%    d/dt psi_f = r_d + K_f * ssa( psi_d - psi_f ) )
%
% where the desired yaw rate r_d = d(psi_d)/dt is computed by
%
%    d/dt y_e = -U * y_e / sqrt( Delta^2 + y_e^2 )
%    r_d = -Kp * dy_e/dt / ( 1 + (Kp * y_e)^2 )   
%
% Initialization:
%   The active waypoint (xk, yk) where k = 1,2,...,n is a persistent
%   integer should be initialized to the first waypoint, k = 1, using 
%   >> clear ALOSpsi
%
% Inputs:   
%   (x,y): craft North-East positions (m)
%   Delta: positive look-ahead distance (m)
%   gamma: positive adaptive gain constant
%   h: sampling time (s)
%   R_switch: go to next waypoint when the along-track distance x_e 
%             is less than R_switch (m)
%   wpt.pos.x = [x1, x2,..,xn]' array of waypoints expressed in NED (m)
%   wpt.pos.y = [y1, y2,..,yn]' array of waypoints expressed in NED (m)
%   U: speed, vehicle cruise speed or time-varying measurement (m/s)
%   K_f: observer gain for desired yaw angle (typically 0.1-0-5)
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%      dist = sqrt(  (wpt.pos.x(k+1)-wpt.pos.x(k))^2 
%                  + (wpt.pos.y(k+1)- wpt.pos.y(k))^2 );
%
% Outputs:  
%    psi_d: desired heading angle (rad)
%    r_d:   desired yaw rate (rad/s)
%    y_e:   cross-track error (m)
%
% For course control use the functions LOSchi.m and ILOSchi.m.
%
% Ref. Fossen, T. I. (2023). An Adaptive Line-of-sight (ALOS) Guidance Law 
% for Path Following of Aircraft and Marine Craft. IEEE Transactions on 
% Control Systems Technology, 1-8. 
% https://doi.org/10.1109/TCST.2023.3259819
%  
% Author:    Thor I. Fossen
% Date:      23 April 2023
% Revisions: 

persistent k;        % active waypoint index (initialized by: clear ALOSpsi)
persistent xk yk;    % active waypoint (xk, yk) corresponding to integer k
persistent psi_f;    % filtered heading angle command
persistent beta_hat; % estimate of the crab angle

%% Initialization of (xk, yk) and (xk_next, yk_next), and integral state 
if isempty(k)

    % check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end

    % check input parameters
    if (R_switch < 0); error("R_switch must be larger than zero"); end
    if (Delta < 0); error("Delta must be larger than zero"); end

    beta_hat = 0;           % initial states
    psi_f = 0;
    k = 1;                  % set first waypoint as the active waypoint
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
    bearing = atan2((wpt.pos.y(n)- wpt.pos.y(n-1)), (wpt.pos.x(n)-wpt.pos.x(n-1)));
    R = 1e10;
    xk_next = wpt.pos.x(n) + R * cos(bearing);
    yk_next = wpt.pos.y(n) + R * sin(bearing); 
end

%% Compute the path-tangnetial angle w.r.t. North
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

% ALOS guidance law
Kp = 1/Delta;
psi_d = psi_f;
psi_ref = pi_h - beta_hat - atan( Kp * y_e ); 

% desired yaw rate r_d
Dy_e = -U * y_e / sqrt( Delta^2 + y_e^2 );    % d/dt y_e
r_d = -Kp * Dy_e / ( 1 + (Kp * y_e)^2 );      % d/dt psi_d

% integral state: y_int[k+1]
beta_hat = beta_hat + h * gamma * Delta * y_e / sqrt( Delta^2 + y_e^2 );

% observer/filter: psi_f[k+1]
psi_f = psi_f + h * ( r_d + K_f * ssa( psi_ref - psi_f ) );


end

