function SIMfrigate()
% SIMfrigate is compatibel with MATLAB and GNU Octave (www.octave.org)
%
% User-editable script for simulation of a heading controlled frigate
% (Length 100 m) using PID control with reference feedforward. The frigate 
% is modelled by the Norrbin (1963) nonlinear model.
%
% Calls:      
%   frigate.m   : Frigate dynamics
%   refModel.m  : Reference model dynamics
%
% Author:     Thor I. Fossen
% Date:       2024-04-20
% Revisions:  

close all;

%% USER INPUTS
h  = 0.05;                      % sample time (s)
N  = 10000;                     % number of samples

%% AUTOPILOTS PARAMETERS
U = 10;                         % speed (m/s)
delta_max = deg2rad(30);        % max rudder angle (deg)

% Reference model parameters
wn_d = 0.1;                     % natural frequency (rad/s)
zeta_d = 1.0;                   % relative damping factor (-)
r_max = deg2rad(5.0);           % maximum turning rate (rad/s)

% PID heading autopilot 
T = 27.0;                       % Nomoto gains at U = 9 m/s
K = 0.18;
wn = 0.5;                       % closed-loop natural frequency (rad/s)
zeta = 1.0;                     % closed loop relative damping factor (-)

Kp = (T/K) * wn^2;                    
Kd = (T/K) * (2 * zeta * wn - 1/T);
Td = Kd / Kp; 
Ti = 10 / wn;

% Initial states
psi = 0;                        % heading angle (rad)
r = 0;                          % yaw rate (rad/s)
delta = 0;                      % rudder angle (rad)
psi_int = 0;                    % integral state
psi_d = psi;                    % reference model states
r_d = r;
a_d = 0;

%% MAIN LOOP
simdata = zeros(N+1,7);              % memory allocation

for i = 1:N+1

    t= (i-1) * h;                    % simulation time in seconds

    % Reference model, step input
    psi_ref = deg2rad(30);
    if t > 100
        psi_ref = deg2rad(70);
    end

   % Reference model propagation
    [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
            zeta_d, wn_d, h, 1);

    delta_c = (T/K) * a_d + (1/K) * r_d -...
        Kp * (ssa( psi-psi_d) +...
        Td * (r - r_d) + (1/Ti) * psi_int );

    delta_c = sat(delta_c, delta_max);   % amplitude saturation

    % Frigate dynamics
    [psi_dot, r_dot, delta_dot] = frigate(r, U, delta, delta_c);

    % Store data for presentation
    simdata(i,:) = [t, psi, r, delta, delta_c, psi_d, r_d];

    % Euler's integration method (k+1)
    psi = psi + h * psi_dot;
    r = r + h * r_dot;
    delta = delta + h * delta_dot;
    psi_int = psi_int + h * ssa( psi - psi_d );

end

%% PLOTS
% simdata(i,:) = [t, psi, r, delta, delta_c, psi_d, r_d]
t       = simdata(:,1);      
psi     = rad2deg(simdata(:,2)); 
r       = rad2deg(simdata(:,3));          
delta   = rad2deg(simdata(:,4));   
delta_c = rad2deg(simdata(:,5));
psi_d   = rad2deg(simdata(:,6));
r_d   = rad2deg(simdata(:,7));

figure(1)

subplot(221)
plot(t,r,t,r_d)
xlabel('Time (s)')
title('Yaw rate (deg/s)')
legend('True','Desired')
grid

subplot(222)
plot(t,psi,t,psi_d)
xlabel('Time (s)')
title('Yaw angle (deg)')
legend('True','Desired')
grid

subplot(212)
plot(t,delta,t,delta_c)
xlabel('Time (s)')
title('Rudder angle \delta (deg)')
legend('True','Commanded')
grid

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

end

