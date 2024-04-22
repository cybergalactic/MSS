function [psi_ref,y_e] = ALOSpsi(x,y,Delta_h,gamma_h,h,R_switch,wpt)
% ALOSpsi is compatible with MATLAB and GNU Octave (www.octave.org). The
% function [psi_ref,y_e] = ALOSpsi(x,y,Delta_h,gamma_h,h,R_switch,wpt) 
% ALOSpsi computes the desired heading angle, psi_ref, and cross-track
% error, y_e, for paths that consist of straight lines connecting waypoints 
% (wpt.pos.x, wpt.pos.y). The desired heading angle is computed using the 
% adaptive LOS (ALOS) guidance law by Fossen (2023),
%
%  psi_ref[k] = pi_h - beta_hat[k] - atan( y_e[k] / Delta_h )
%
%  beta_hat[k+1] = beta_hat[k] + h * gamma_h * Delta_h * ...
%       y_e[k] / sqrt( Delta_h^2 + y_e[k]^2 )
%
% where pi_h is the path-tangential (azimuth) angle with respect to the 
% North axis and y_e is the cross-track error. To handle steps in the yaw 
% angle command, psi_ref, due to waypoint switching, the following observer 
% should be used:
%
%  [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f, T_f)
%
% where psi_d is the desired yaw angle and r_d is the desired yaw rate r_d. 
%
% Initialization:
%   The active waypoint (xk, yk) where k = 1,2,...,n is a persistent
%   integer should be initialized to the first waypoint, k = 1, using 
%   >> clear ALOSpsi
%
% Inputs:   
%   (x,y): craft North-East positions (m)
%   Delta_h: positive look-ahead distance (m), (typically 5-20 m)
%   gamma_h: positive adaptive gain constant, (typically 0.001)
%   h: sampling time (s)
%   R_switch: go to next waypoint when the along-track distance x_e 
%             is less than R_switch (m)
%   wpt.pos.x = [x1, x2,..,xn]' array of waypoints expressed in NED (m)
%   wpt.pos.y = [y1, y2,..,yn]' array of waypoints expressed in NED (m)
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%     dist = sqrt(   ( wpt.pos.x(k+1) - wpt.pos.x(k) )^2 
%                  + ( wpt.pos.y(k+1) - wpt.pos.y(k) )^2 )
% Outputs:  
%    psi_ref: LOS heading angle (rad)
%    y_e:     cross-track error (m)
%
% Reference:
%   T. I. Fossen (2023). An Adaptive Line-of-sight (ALOS) Guidance Law 
%   for Path Following of Aircraft and Marine Craft. IEEE Transactions on 
%   Control Systems Technology 31(6),2887-2894. 
%   https://doi.org/10.1109/TCST.2023.3259819
%
% For course control use the functions LOSchi.m and ILOSchi.m.
% See also: SIMotter.m and demoOtterUSVPathFollowingHeadingControl.slx
% for example implementations using an Otter USV.
%  
% Author:    Thor I. Fossen
% Date:      2023-03-23
% Revisions: 2024-04-01 - Added compability to LOSobserver.m

persistent k;        % active waypoint index, initialized by 'clear ALOSpsi'
persistent xk yk;    % active waypoint (xk, yk) corresponding to integer k
persistent beta_hat; % estimate of the crab angle

%% Initialization of (xk, yk) and (xk_next, yk_next), and integral state 
if isempty(k)

    % Check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end

    % Check input parameters
    if (R_switch < 0); error("R_switch must be larger than zero"); end
    if (Delta_h < 0); error("Delta must be larger than zero"); end

    beta_hat = 0;           % initial states
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

% Along-track and cross-track errors (x_e, y_e) 
x_e =  (x-xk) * cos(pi_h) + (y-yk) * sin(pi_h);
y_e = -(x-xk) * sin(pi_h) + (y-yk) * cos(pi_h);

% If the next waypoint satisfy the switching criterion, k = k + 1
d = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
if ( (d - x_e < R_switch) && (k < n) )
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
    fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);
end

% ALOS guidance law
psi_ref = pi_h - beta_hat - atan( y_e/Delta_h ); 

% Propagation of states to time k+1
beta_hat = beta_hat + h * gamma_h * Delta_h * y_e / sqrt( Delta_h^2 + y_e^2 );

end

