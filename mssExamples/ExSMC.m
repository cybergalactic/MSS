% ExSMC integral sliding mode control (SMC) design for heading autopilot. 
% A conventional SMC is designed for the Norrbin (1963) nonlinear yaw model
% 
%                       psi_dot = r
%       T r_dot + n3 r^3 + n1 r = K delta + d_r
%
% where d_r is an unknown bounded disturbance in yaw and 
%
%   [psi_dot, r_dot, delta_dot] = ROVzefakkel(r,U,delta,delta_c,d_r)
%
% is the Norrbin model for the ROV Zefakkel (Length 45 m) inclduing
% actuator dynamics and saturation
%
% Author:    Thor I. Fossen
% Date:      19 June 20200
% Revisions: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h    = 0.05;     % sampling time [s]
N  =   6000;    % no. of samples

psi_ref = 10 * pi/180;  % desired yaw angle

flag = 3;      % 1 = conventional SMC using sgn(sigma)
               % 2 = conventional SMC using tanh(sigma/phi)
               % 3 = conventional SMC using sat(sigma)
               
% SMC parameters
T_hat = 30;             % ship time constant
K_hat = 0.3;            % ship gain constant 
n3 = 0.4;               % 3rd-order maneuvering coefficient
n1 = 1;                 % 1st-order maneuvering coeffiecient, stable ship
K_sigma = 0.1;          % SMC gain
phi = 0.001;            % boundary layer parameter tanh(sigma/phi)
lambda = 0.1;           % sliding variable parameter lambda > 0
z_psi = 0;              % intial integral state

% ship model parameters
psi = 0;                % initial yaw angle (rad)
r = 0;                  % initial yaw rate (rad/s)
delta = 0;              % initial rudder angle (rad)
U = 4;                  % ship cruise speed (m/s)
d_r = 0.37*(1*pi/180);  % d_r = K delta_0 (1 degree unknown rudder bias)

% reference model
wn = 0.1;               % reference model nataural frequnecy
xd = [ 0 0 0]';         % inital reference model states (3rd order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(N+1,8);                  % table of simulation data

for i=1:N+1

    t = (i-1) * h;                      % time (s)   
    
    if (i > 2000), psi_ref = -10 * pi/180;  end 
    if (i > 2000 && i>4000), psi_ref = 20 * pi/180;  end 
    
    % 3rd-order reference model for yaw
    Ad = [ 0 1 0
           0 0 1 
           -wn^3 -3*wn^2 -3*wn ];
    Bd = [0 0 wn^3]';
           
    xd_dot = Ad * xd + Bd * psi_ref;    
    
    % conventional sliding mode controller 
    % the rudder rudder dynamics is unknown and |d_r | < d_r^max
    e_psi   = ssa( psi - xd(1) );
    e_r     = r - xd(2);
    sigma = e_r + 2 * lambda * e_psi + lambda^2 * z_psi; 
    r_r = r - sigma;  % sigma = r - r_r
    r_r_dot = xd_dot(3) - 2 * lambda * e_r - lambda^2 * e_psi;
    
    switch flag
        case flag==1        % sgm(sigma)
            delta_c = (1/K_hat) * ( T_hat * r_r_dot + (n3 * r^2 + n1) * r_r... 
            - K_sigma * sign(sigma) );
        case flag == 2      % tanh(sigma)
            delta_c = (1/K_hat) * ( T_hat * r_r_dot + (n3 * r^2 + n1) * r_r... 
            - K_sigma * tanh(sigma/phi) );
        otherwise           % sat(sigma)
            if (abs(sigma/phi) > 1)
              delta_c = (1/K_hat) * ( T_hat * r_r_dot + (n3 * r^2 + n1) * r_r... 
              - K_sigma * sign(sigma/phi) );
            else
              delta_c = (1/K_hat) * ( T_hat * r_r_dot + (n3 * r^2 + n1) * r_r... 
              - K_sigma * (sigma/phi) );
            end
    end
    
    % Norrbin model for the ROV Zefakkel
    [psi_dot, r_dot, delta_dot] = ROVzefakkel(r,U,delta,delta_c,d_r); 
    % store simulation data in a table (for testing)
    simdata(i,:) = [t psi r delta delta_c xd'];       
     
    % Euler integration
    xd = xd + h * xd_dot;
    psi = psi + h * psi_dot;
    r = r + h * r_dot; 
    delta = delta + h * delta_dot; 
    z_psi = z_psi + h * e_psi;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);           
psi     = (180/pi) * simdata(:,2); 
r       = (180/pi) * simdata(:,3); 
delta   = (180/pi) * simdata(:,4);
delta_c = (180/pi) * simdata(:,5);           
psi_d   = (180/pi) * simdata(:,6); 
r_d     = (180/pi) * simdata(:,7);

figure(gcf)
subplot(311)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(312)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');
