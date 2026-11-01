function x_prd = EKF_5states(position1, position2,...
    h, Z, frame, Qd, Rd, alpha_1, alpha_2, x_prd_init) 
% This function estimates the Speed Over Ground (SOG), Course Over Ground 
% (COG), and course rate from GNSS positions measurements (xn[k], yn[k]) 
% expressed in NED or latitude-longitude (mu[k], l[k]) using the 5-state 
% discrete-time extended Kalman filter (EKF). The output is the next-step
% prior estimate x_prd[k+1], including positions (x, y), SOG (U_g), 
% COG (chi), and course rate (omega_chi). Call at the end of the loop with
% measurements from sample k, then use the result at sample k+1.
% The EKF discrete-time state-space 
% model is (Fossen and Fossen, 2021):
%
%   x[k+1] = x[k] + h * U[k] * cos(chi[k])
%   y[k+1] = y[k] + h * U[k] * sin(chi[k])
%   U[k+1] = U[k] + h * ( -alpha_1 * U[k] + w_1[k] )  
%   chi[k+1] = chi[k] + h * omega_chi[k]
%   omega_chi[k+1] = omega_chi[k] + h * ( -alpha_2 * omega_chi[k] + w_2[k] ) 
%
%   y_1[k] = x[k] + v_1[k]
%   y_2[k] = y[k] + v_2[k]
% 
% The script calling EKF_5states.m should include the command: 
% 
%   clear LOSchi EKF_5states  
% 
% to clear the persistent variables x_prd_persistent, P_prd, and count.
%
% Examples:
%   x_prd = EKF_5states(x, y, h, Z, 'NED', Qd, Rd) 
%   x_prd = EKF_5states(x, y, h, Z, 'NED', Qd, Rd, alpha_1, alpha_2, x_prd_init) 
% where x_prd = [x, y, U, chi, omega_chi] 
%
%   x_prd = EKF_5states(mu, l, h, Z, 'LL', Qd, Rd)
%   x_prd = EKF_5states(mu, l, h, Z, 'LL', Qd, Rd, alpha_1, alpha_2, x_prd_init)
% where x_prd = [mu, l, U, chi, omega_chi]
%
% Inputs:
%  position1, position2: North-East positions (m) or Latitude-Longitude (rad)
%  h:         Sampling time (s)
%  Z:         Position measurement frequency (Z times slower). Must be integer
%  frame:     'NED' (North-East-Down) or 'LL' (Latitude-Longitude)
%  Qd:        2x2 covariance of sampled process noise w[k]; Ed includes h
%  Rd:        EKF 2x2 position measurement covariance matrix
%  alpha_1:   (Optionally), Singer constant, speed
%  alpha_2:   (Optionally), Singer constant, course rate
%  x_prd_init: (Optionally), Initial predicted state vector x_prd_persistent
%
% See exOtterEKF.m and SIMmariner.m for case studies using EKF_5states.m to 
% estimate the COG, SOG, and course rate.
%
% Simulink Models:
%   demoMarinerPathFollowingCourseControl.slx : Simulink model demonstrating  
%       path-following course control using straight lines and waypoint 
%       switching.
%
% References: 
%   S. Fossen and T. I. Fossen (2021). Five-state Extended Kalman Filter for 
%     Estimation of Speed Over Ground (SOG), Course Over Ground (COG) and
%     Course Rate of Unmanned Surface Vehicles (USVs): Experimental Results. 
%     Sensors 21(23). 
%
% Author:   Thor I. Fossen
% Date:     2021-07-25 
% Revisions:
%   2024-07-27 - Added optional initial state vector x_prd_init.
%   2026-11-01 - Include WGS-84 curvature-radius derivatives in the LL
%                Jacobian and allow separate optional Singer constants.
%                Return the next-step prediction, consistently with MSS INS
%                functions, and clarify sampled-noise covariance.

persistent x_prd_persistent;
persistent P_prd;
persistent count;

I5 = eye(5);

if isempty(x_prd_persistent) 
    if nargin < 10  % Check if x_prd_init is provided
        disp(['Using default initial EKF states: x_prd = ' ...
            '[position1 position2 0 0 0]']);
        x_prd_persistent = [position1 position2 0 0 0]';
    else
        x_prd_init = x_prd_init(:);
        disp(['Using user specified initial EKF states: x_prd = [', ...
            num2str(x_prd_init', '%.2f '), ']']);
       x_prd_persistent = x_prd_init;
    end
       P_prd = I5;          
       count = 1;
end

% WGS-84 data
a = 6378137;           % Semi-major axis (equatorial radius)
f = 1/298.257223563;   % Flattening 
e = sqrt(2*f - f^2);   % Earth eccentricity

if nargin < 8
    alpha_1 = 0.01;        % Default Singer constant, speed
end
if nargin < 9
    alpha_2 = 0.1;         % Default Singer constant, course rate
end

Cd  = [1 0 0 0 0            
       0 1 0 0 0 ];                                    
Ed = h * [ 0 0; 0 0; 1 0; 0 0; 0 1 ]; 

if (count == 1)   
    % Correct estimates if new measurement
    y  = [ position1        % y1 = x^n    or   y1 = latitude
           position2 ];     % y2 = y^n    or   y2 = longitude
              
    K = P_prd * Cd' / ( Cd * P_prd * Cd' + Rd );   % KF gain 
    IKC = I5 - K * Cd;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';      % corrector (k)
    eps = y - Cd * x_prd_persistent; 
    if strcmp(frame,'LL'); eps = ssa(eps); end
    x_hat = x_prd_persistent + K * eps; 
    count = Z;            
else                    
    % Predict/propagate states and covariance if no new measurements
    x_hat = x_prd_persistent;
    P_hat = P_prd; 
    count = count - 1;   
end

if strcmp(frame,'NED')  % x = [ x^n y^n U chi omega_chi ]'
   
    f = [ x_hat(3) * cos(x_hat(4))  
          x_hat(3) * sin(x_hat(4))
         -alpha_1 * x_hat(3)               
          x_hat(5)
         -alpha_2 * x_hat(5) ];
  
    Ad = I5 + h * ...
        [ 0 0 cos(x_hat(4)) -x_hat(3)*sin(x_hat(4)) 0 
          0 0 sin(x_hat(4))  x_hat(3)*cos(x_hat(4)) 0 
          0 0 -alpha_1 0 0
          0 0 0 0 1
          0 0 0 0 -alpha_2];
    
elseif strcmp(frame,'LL')  % x = [ mu l U chi omega_chi ]'
    
    Rn = a / sqrt( 1-e^2 * sin(x_hat(1))^2 );
    Rm = Rn * ( (1-e^2) / (1-e^2 * sin(x_hat(1))^2) );
    % Logarithmic derivatives of the curvature radii with respect to latitude.
    dlogRn = e^2 * sin(x_hat(1)) * cos(x_hat(1)) ...
        / (1-e^2 * sin(x_hat(1))^2);
    dlogRm = 3 * dlogRn;
    
    f = [ ( 1 / Rm ) * x_hat(3) * cos(x_hat(4))  
          ( 1 / ( Rn * cos(x_hat(1)) ) ) * x_hat(3) * sin(x_hat(4))
         -alpha_1 * x_hat(3)               
          x_hat(5)
         -alpha_2 * x_hat(5) ];
    
    Ad = I5 + h * ...
        [ -dlogRm/Rm*x_hat(3)*cos(x_hat(4)) ...
                0 (1/Rm)*cos(x_hat(4)) -(1/Rm)*x_hat(3)*sin(x_hat(4)) 0 
          (tan(x_hat(1))-dlogRn) / (Rn * cos(x_hat(1))) ...
                * x_hat(3) * sin(x_hat(4)) ...
                0 (1/(Rn*cos(x_hat(1))))*sin(x_hat(4))...
               (1/(Rn*cos(x_hat(1)))) * x_hat(3) * cos(x_hat(4)) 0 
          0 0 -alpha_1 0 0
          0 0 0 0 1
          0 0 0 0 -alpha_2];    
end   

% Predictor (k+1)   
x_prd = x_hat + h * f;
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

% Return the prior prediction for the next sampling instant.
x_prd_persistent = x_prd;

end
