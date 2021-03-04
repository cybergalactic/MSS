% ExSTA adaptive-gain super twisting algorithm (STA) for heading control. 
% The yaw dynamics is based on the Norrbin (1963) nonlinear model
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
% Date:      20 June 2020
% Revisions: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.05;    % sampling time [s]
N  = 8000;    % no. of samples

psi_ref = 10 * pi/180;  % desired yaw angle
               
% STW parameters
T = 30;                 % ship time constant
K = 0.3;                % ship gain constant
lambda = 0.1;           % sliding variable parameter lambda > 0
alpha_0 = 0.03;         % adaptation gain
beta_0 = 0.0001;       % adaptation gain
v = 0;                  % initial state
alpha = 0;              % initial state

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
simdata = zeros(N+1,11);                % table of simulation data

for i=1:N+1

    t = (i-1) * h;                      % time (s)   
    
    if (i > 2000), psi_ref = -10 * pi/180;  end 
    if (i > 4000), psi_ref =  10 * pi/180;  end 
    if (i > 6000), psi_ref = -10 * pi/180;  end    
    
    % 3rd-order reference model for yaw
    Ad = [ 0 1 0
           0 0 1 
           -wn^3 -3*wn^2 -3*wn ];
    Bd = [0 0 wn^3]';
           
    xd_dot = Ad * xd + Bd * psi_ref;    
    
    % sliding variable with integral action 
    e_psi   = ssa( psi - xd(1) );
    e_r     = r - xd(2);
    sigma = e_r + lambda * e_psi; 
    
    % STW sliding mode controller    
    if (abs(sigma) < 0.01)
        alpha_dot = 0;
    else
        alpha_dot = alpha_0;
    end
    
    beta = beta_0;
    phi = 0.01;
    v_dot = -beta * tanh(sigma/phi);  
    w = -alpha * sqrt(abs(sigma)) * sign(sigma) + v;
    delta_c = (T/K) * w;
    
    % Norrbin model for the ROV Zefakkel
    [psi_dot,r_dot,delta_dot] = ROVzefakkel(r,U,delta,delta_c,d_r); 
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t psi r delta delta_c xd' alpha beta v];       
     
    % Euler integration
    xd = xd + h * xd_dot;
    psi = psi + h * psi_dot;
    r = r + h * r_dot; 
    delta = delta + h * delta_dot; 
    v = v + h * v_dot;
    alpha = alpha + h * alpha_dot;
    
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
a_d     = (180/pi) * simdata(:,8);
alpha   = simdata(:,9);
beta    = 100 * simdata(:,10);
v       = 100 * simdata(:,11);

figure(gcf)
subplot(411)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(412)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');
subplot(413)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');
subplot(414)
plot(t,alpha,'r',t,beta,'b:',t,v,'k-.','linewidth',2);
title('STW variables'); xlabel('time (s)');
legend('\alpha','100 \beta','100 v')

