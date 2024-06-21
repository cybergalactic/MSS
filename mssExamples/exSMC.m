% ExSMC is compatible with MATLAB and GNU Octave (www.octave.org).
% Integral sliding mode control (SMC) design for heading autopilot. 
% A conventional SMC is designed for the Norrbin (1963) nonlinear yaw model
% 
%                       psi_dot = r
%       T r_dot + n3 r^3 + n1 r = K delta + d_r
%
% where d_r is an unknown bounded disturbance in yaw and 
%
%   [psi_dot, r_dot, delta_dot] = zeefakkel(r,U,delta,delta_c,d_r)
%
% is the Norrbin model for the Zeefakkel (Length 45 m), including
% actuator dynamics and saturation.
%
% Author:    Thor I. Fossen
% Date:      19 June 20200
% Revisions: 30 Aug 2023 - minor updates and improvements

%% USER INPUTS
h  = 0.05;    % sampling time [s]
N  = 6000;    % no. of samples

psi_ref = deg2rad(10);  % desired yaw angle

flag = 3;     % 1 = conventional SMC using sgn(sigma)
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
d_r = 0.37*deg2rad(1);  % d_r = K delta_0 (1 degree unknown rudder bias)

% reference model
wn = 0.1;               % reference model nataural frequnecy
xd = [0 0 0]';          % inital reference model states (3rd order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(N+1,8);  % table of simulation data

for i=1:N+1

    t = (i-1) * h;       % time (s)   
    
    if (t > 100), psi_ref = deg2rad(-10); end 
    if (t > 200), psi_ref = deg2rad(20);  end 
    
    % 3rd-order reference model for yaw
    Ad = [ 0 1 0
           0 0 1 
           -wn^3 -3*wn^2 -3*wn ];
    Bd = [0 0 wn^3]';
           
    xd_dot = Ad * xd + Bd * psi_ref;    
    
    % Conventional sliding mode controller 
    % The rudder rudder dynamics is unknown and |d_r | < d_r^max
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
    [psi_dot, r_dot, delta_dot] = zeefakkel(r,U,delta,delta_c,d_r); 
    
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
psi     = rad2deg(simdata(:,2)); 
r       = rad2deg(simdata(:,3)); 
delta   = rad2deg(simdata(:,4));
delta_c = rad2deg(simdata(:,5));           
psi_d   = rad2deg(simdata(:,6)); 
r_d     = rad2deg(simdata(:,7));

figure(gcf)
subplot(311)
plot(t,psi,t,psi_d);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(312)
plot(t,r,t,r_d);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
