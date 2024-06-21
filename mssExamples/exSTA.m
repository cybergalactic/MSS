% ExSTA is compatible with MATLAB and GNU Octave (www.octave.org).
% Super twisting adaptive-gain (STA) sliding mode control for heading
% control. The yaw dynamics is based on the Norrbin (1963) nonlinear model
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
% Date:      20 June 2020
% Revisions: 29 Nov 2022 - Beta is changed to beta = beta_1 * alpha + beta_0
%            30 Aug 2023 - Minor updates and improvements

%% USER INPUTS
h  = 0.01;     % sampling time (s)
N  = 40000;    % no. of samples

psi_ref = deg2rad(10);  % desired yaw angle
               
% STW parameters
T = 30;                 % ship time constant
K = 0.3;                % ship gain constant
lambda = 0.05;          % sliding variable parameter lambda > 0
alpha_0 = 0.01;         % adaptation gains
beta_0 = 0.00001;       
beta_1 = 0.0001;

v = 0;                  % initial STA states
alpha = 0;              

% ship model parameters
psi = 0;                % initial yaw angle (rad)
r = 0;                  % initial yaw rate (rad/s)
delta = 0;              % initial rudder angle (rad)
U = 4;                  % ship cruise speed (m/s)
d_r = K * deg2rad(1);   % d_r = K delta_0 (1 degree unknown rudder bias)

% reference model
wn = 0.1;               % reference model nataural frequnecy
xd = [0 0 0]';          % inital reference model states (3rd order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(N+1,11);                % table of simulation data

for i=1:N+1

    t = (i-1) * h;                      % time (s)   
    
    % setpoints and steps in rudder bias 
    if (t > 100), psi_ref = -10 * pi/180; d_r =  2*K * deg2rad(1); end 
    if (t > 200), psi_ref =  10 * pi/180; d_r = -2*K * deg2rad(1); end 
    if (t > 300), psi_ref = -10 * pi/180; d_r = 0; end    
    
    % 3rd-order reference model for yaw
    Ad = [ 0 1 0
           0 0 1 
           -wn^3 -3*wn^2 -3*wn ];
    Bd = [0 0 wn^3]';
           
    xd_dot = Ad * xd + Bd * psi_ref;    
    
    % sliding variable 
    e_psi   = ssa( psi - xd(1) );
    e_r     = r - xd(2);
    sigma = e_r + lambda * e_psi; 
    
    % STW adaptive-gain sliding mode controller    
    if (abs(sigma) < 0.001)
        alpha_dot = 0;
    else
        alpha_dot = alpha_0;
    end
    
    beta = beta_1 * alpha + beta_0;
    v_dot = -beta * sign(sigma); 
    w = -alpha * sqrt(abs(sigma)) * sign(sigma) + v;
    delta_c = (T/K) * w;
    
    % Norrbin model for the Zefakkel
    [psi_dot,r_dot,delta_dot] = zeefakkel(r,U,delta,delta_c,d_r); 
    
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
psi     = rad2deg(simdata(:,2)); 
r       = rad2deg(simdata(:,3)); 
delta   = rad2deg(simdata(:,4));
delta_c = rad2deg(simdata(:,5));           
psi_d   = rad2deg(simdata(:,6)); 
r_d     = rad2deg(simdata(:,7));
a_d     = rad2deg(simdata(:,8));
alpha   = simdata(:,9);
beta    = 1000 * simdata(:,10);
v       = 100 * simdata(:,11);

figure(gcf)
subplot(411)
plot(t,psi,t,psi_d);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(412)
plot(t,r,t,r_d);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');
subplot(413)
plot(t,delta,t,delta_c);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');
subplot(414)
plot(t,alpha,'r',t,beta,'b:',t,v,'k-.');
title('STW variables'); xlabel('time (s)');
legend('\alpha','1000 \beta','100 v')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

