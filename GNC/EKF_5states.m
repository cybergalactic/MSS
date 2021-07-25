function x_hat = EKF_5states(GNSS1,GNSS2,h_samp,Z,frame,Qd,Rd)  
% EKF_5states estimates SOG, COG, and course rate from GNSS positions
% (xn[k], yn[k]) expressed in NED or latitude-longitude (mu[k], l[k]) using 
% a 5-states discrete-time extended Kalamn filter (EKF). The output is the 
% predicted state vector x_hat[k+1] as defined below.
%
% Outputs:
%  [x,y,U,chi omega_chi]  = EKF_5states(x,y,h,Z,'NED',Qd,Rd) 
%  [mu,l,U,chi,omega_chi] = EKF_5states(mu,l,h,Z,'LL',Qd,Rd)
%
% Inputs:
%  GNSS1,GNSS2: North-East positions (m) or Latitude-Longitude (rad)
%  h_samp:      EKF sampling time (s), h_samp = 1 / measurement frequency
%  Z:           h_samp * GNSS measurement frequency (Hz) (must be integer)
%  frame:       'NED' (North-East) or 'LL' (Latitude-Longitude)
%  Qd:          EKF 2x2 process cov. matrix for speed and course rate
%  Rd:          EKF 2x2 GNSS measurement cov. matrix
%
% Reference: S. Fossen and T. I. Fossen (2021). Five-state extended Kalman 
% filter for estimation of speed over ground(SOG), course over ground (COG) 
% and course rate of surface vehicles. Journal of Marine Science and 
% Applications. Submitted.
%
% Author:   Thor I. Fossen
% Date:     25 July 2021
% Revisions: 

persistent x_prd;
persistent P_prd;
persistent count;

I5 = eye(5);

if isempty(x_prd)
   disp('Init EKF states') 
   x_prd = [ GNSS1 GNSS2 0 0 0 ]';
   P_prd = I5;          
   count = 1;
end

% WGS-84 data
a = 6378137;           % Semi-major axis (equatorial radius)
f = 1/298.257223563;   % Flattening 
e = sqrt(2*f - f^2);   % Earth eccentricity

alpha_1 = 0.01;        % Singer constant, speed
alpha_2 = 0.1;         % Singer constant, course rate

Cd  = [1 0 0 0 0            
       0 1 0 0 0 ];                                    
Ed = h_samp * [ 0 0; 0 0; 1 0; 0 0; 0 1 ]; 

if (count == 1)         % update if new measurement
    y  = [ GNSS1        % y1 = x^n    or   y1 = latitude
           GNSS2 ];     % y2 = y^n    or   y2 = longitude
              
    K = P_prd * Cd' / ( Cd * P_prd * Cd' + Rd );   % KF gain 
    IKC = I5 - K * Cd;
    P_hat = IKC * P_prd * IKC' + K * Rd * K';      % corrector (k)
    eps = y - Cd * x_prd; 
    if strcmp(frame,'LL'); eps = ssa(eps); end
    x_hat = x_prd + K * eps; 
    count = Z;            
else                    % no update
    x_hat = x_prd;
    P_hat = P_prd; 
    count = count - 1;   
end

if strcmp(frame,'NED')  % x = [ x^n y^n U chi omega_chi ]'
   
    f = [ x_hat(3) * cos(x_hat(4))  
          x_hat(3) * sin(x_hat(4))
         -alpha_1 * x_hat(3)               
          x_hat(5)
         -alpha_2 * x_hat(5) ];
  
    Ad = I5 + h_samp * ...
        [ 0 0 cos(x_hat(4)) -x_hat(3)*sin(x_hat(4)) 0 
          0 0 sin(x_hat(4))  x_hat(3)*cos(x_hat(4)) 0 
          0 0 -alpha_1 0 0
          0 0 0 0 1
          0 0 0 0 -alpha_2];
    
elseif strcmp(frame,'LL')  % x = [ mu l U chi omega_chi ]'
    
    Rn = a / sqrt( 1-e^2 * sin(x_hat(1))^2 );
    Rm = Rn * ( (1-e^2) / (1-e^2 * sin(x_hat(1))^2) );
    
    f = [ ( 1 / Rm ) * x_hat(3) * cos(x_hat(4))  
          ( 1 / ( Rn * cos(x_hat(1)) ) ) * x_hat(3) * sin(x_hat(4))
         -alpha_1 * x_hat(3)               
          x_hat(5)
         -alpha_2 * x_hat(5) ];
    
    Ad = I5 + h_samp * ...
        [ 0 0 (1/Rm)*cos(x_hat(4)) -(1/Rm)*x_hat(3)*sin(x_hat(4)) 0 
          tan(x_hat(1)) / (Rn * cos(x_hat(1))) * x_hat(3) * sin(x_hat(4))...
                0 (1/(Rn*cos(x_hat(1))))*sin(x_hat(4))...
               (1/(Rn*cos(x_hat(1)))) * x_hat(3) * cos(x_hat(4)) 0 
          0 0 -alpha_1 0 0
          0 0 0 0 1
          0 0 0 0 -alpha_2];    
end   

% Predictor (k+1)   
x_prd = x_hat + h_samp * f;
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

end
