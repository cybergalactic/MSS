function [psi_ref, theta_ref, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
    ALOS3D(x,y,z,Delta_h,Delta_v,gamma_h,gamma_v,M_theta,h,R_switch,wpt)
% ALOS3D is compatible with MATLAB and GNU Octave (www.octave.org). The
% function [psi_ref, theta_ref, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
% ALOS3D(x,y,z,Delta_h,Delta_v,gamma_h,gamma_v,M_theta,h,R_switch,wpt)
% computes the desired heading angle psi_d and pitch angle theta_d when
% the path is a straight line segment going through the waypoints 
% (wpt.pos.x,wpt.pos.y,wpt.pos.z). The desired heading and pitch angles are
% computed using the adaptive LOS (ALOS) guidance law by Fossen and Aguiar
% (2024) where
%
%   theta_ref[k] = pi_v + alpha_hat[k] + atan( z_e[k] / Delta_v )
%   psi_ref[k]   = pi_h - beta_hat[k]  - atan( y_e[k] / Delta_h )
%
%   alpha_hat[k+1] = gamma_v * Delta_v / ...
%       sqrt( Delta_v^2 + z_e[k]^2 ) * Proj( alpha_hat, z_e[k] )
%   beta_hat[k+1]  = gamma_h * Delta_h / ...
%       sqrt( Delta_h^2 + y_e[k]^2 ) * Proj( beta_hat, y_e[k] )
%
% To handle steps in the pitch and yaw angle commands, theta_ref and psi_ref, 
% due to waypoint switching, the following observers should be used:
%
%  [theta_d, q_d] = LOSobserver(theta_d, q_d, theta_ref, h, K_f, T_f)
%  [psi_d, r_d]   = LOSobserver(psi_d, r_d, psi_ref, h, K_f, T_f)
%
% The NED tracking errors are rotated an azimuth angle pi_h about the 
% z-axis and an elevation angle pi_v about the resulting y-axis from the 
% first rotation using
%
%   e = Ry' * Rz' * [x-xk, y-yk, z-zk]'
%
% Initialization:
%   The active waypoint (xk, yk, zk) where k = 1,2,...,n is a persistent
%   integer should be initialized to the first waypoint, k = 1, using 
%   >> clear ALOS3D
%,
% Inputs:   
%   (x,y,z):  craft North-East-Down positions (m)
%   Delta_h:  positive look-ahead distance, horizontal plane (m)
%   Delta_v:  positive look-ahead distance, vertical plane (m)
%   gamma_h:  positive parameter adapation gain, horizontal plane
%   gamma_v:  positive parameter adapation gain, vertical plane
%   M_theta:  maximum parameter: abs(alpha_c) < M_theta, abs(beta_c) < M_theta
%   h:        sampling time (s)
%   R_switch: go to next waypoint when the inside the circle of radius R_switch (m)
%   wpt.pos.x = [x1, x2,..,xn]' array of waypoints expressed in NED (m)
%   wpt.pos.y = [y1, y2,..,yn]' 
%   wpt.pos.z = [z1, z2,..,zn]' 
%
% Feasibility citerion: 
%   The switching parameter R_switch > 0 must satisfy, R_switch < dist, 
%   where dist is the distance between the two waypoints at k and k+1:
%      dist = sqrt(  ( wpt.pos.x(k+1) - wpt.pos.x(k))^2 
%                  + ( wpt.pos.y(k+1) - wpt.pos.y(k))^2 )
%                  + ( wpt.pos.z(k+1) - wpt.pos.z(k))^2 );
% Outputs:  
%    psi_ref:     LOS heading angle (rad)
%    theta_ref:   LOS pitch angle (rad)
%    y_e:         cross-track error (m)
%    z_e:         vertical-track error (m)
%    alpha_c_hat: estimate of alpha_c (rad)
%    beta_c_hat:  estimate of beta_c (rad)
%
% Reference:
% T. I. Fossen and P. Aguiar (2024). A Uniform Semiglobal Exponential 
% Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path Following. 
% Automatica 163, 111556. doi.org/10.1016/j.automatica.2024.111556
%  
% Author:    Thor I. Fossen
% Date:      2022-10-16
% Revisions: 2024-04-01 - Added compability to LOSobserver.m

% Persistent states
persistent k;         % active waypoint index
persistent xk yk zk;  % active waypoint (xk,yk,zk) corresponding to integer k
persistent alpha_hat; % parameter estimate
persistent beta_hat;  % parameter estimate

%% Initialization of (xk,yk,zk) and (xk_next,yk_next,zk_next), and integral state 
if isempty(k)   
  
    % Check if R_switch is smaller than the min. distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ...
            + diff(wpt.pos.z).^2) )
        error("The distances between the waypoints must be larger than R_switch");
    end
    
    % Check input parameters
    if (R_switch < 0)
        error("R_switch must be larger than zero");
    end
    if (Delta_h < 0 || Delta_v < 0)
        error("Delta_h and Delta_v must be larger than zero");
    end         
    
    alpha_hat = 0;        % initial states 
    beta_hat = 0;               

    k = 1;                % set first waypoint as the active waypoint
    xk = wpt.pos.x(k); 
    yk = wpt.pos.y(k);  
    zk = wpt.pos.z(k);     
    fprintf('Active waypoints:\n')
    fprintf('  (x%1.0f,y%1.0f,z%1.0f) = (%.1f,%.1f,%.1f) \n',k,k,k,xk,yk,zk);

end

%% Read next waypoint (xk_next, yk_next, zk_next) from wpt.pos 
n = length(wpt.pos.x);
if k < n                        % if there are more waypoints, read next one 
    xk_next = wpt.pos.x(k+1);  
    yk_next = wpt.pos.y(k+1);   
    zk_next = wpt.pos.z(k+1); 
else                            % else, continue with last bearing and depth
    bearing = atan2((wpt.pos.y(n)-wpt.pos.y(n-1)), (wpt.pos.x(n)-wpt.pos.x(n-1)));
    R = 1e10;
    xk_next = wpt.pos.x(n) + R * cos(bearing);
    yk_next = wpt.pos.y(n) + R * sin(bearing);    
    zk_next = wpt.pos.z(n);
end

%% Compute azimuth ane elevation angles for straigh-line path
pi_h = atan2(  (yk_next-yk) , (xk_next-xk) );  
pi_v = atan2( -(zk_next-zk) , sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 ) );  

% Along-track, cross-track and vertical-track errors (x_e, y_e, z_e) 
Rz = [ cos(pi_h) -sin(pi_h) 0
       sin(pi_h)  cos(pi_h) 0
       0 0 1 ];
Ry = [  cos(pi_v) 0 sin(pi_v)
        0 1 0
       -sin(pi_v) 0 cos(pi_v) ];

e = Ry' * Rz' * [x-xk, y-yk, z-zk]';
y_e = e(2);
z_e = e(3);

% If the next waypoint satisfy the switching criterion, k = k + 1
d = sqrt( (xk_next-x)^2 + (yk_next-y)^2 + (zk_next-z)^2 );
if (d < R_switch) && (k < n)
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
    zk = zk_next;     
    fprintf('  (x%1.0f,y%1.0f,z%1.0f) = (%.1f,%.1f,%.1f) \n',k,k,k,xk,yk,zk);
end

% ALOS guidance laws
psi_ref   = pi_h - beta_hat  - atan( y_e / Delta_h );
theta_ref = pi_v + alpha_hat + atan( z_e / Delta_v );

% ALOS parameter estimates with projection: alpha_hat[k+1], beta_hat[k+1]
alpha_c_hat = alpha_hat;
beta_c_hat  = beta_hat;
alpha_hat = alpha_hat + h * gamma_v * ...
    Delta_v / sqrt( Delta_v^2 + z_e^2 ) * proj(alpha_hat, z_e, M_theta);
beta_hat = beta_hat + h * gamma_h * ...
    Delta_h / sqrt( Delta_h^2 + y_e^2 ) * proj(beta_hat, y_e, M_theta);

end

% *************************************************************************
% Projection algorithm
% *************************************************************************
function y = proj(theta_hat, tau, M_theta)
    
eps = 0.001;
M_theta_hat = M_theta + eps;

if ( abs(theta_hat) > M_theta ) && ( theta_hat * tau > 0 )
    c = min(1, ( M_theta_hat^2 - M_theta^2) / (M_theta_hat^2 - M_theta^2) );
    y = (1 - c) * tau;
else
    y = tau;
end

end


