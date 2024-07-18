function SIMfrigate()
% SIMfrigate is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates a heading-controlled frigate with a length of 
% 100 meters, employing a PID control strategy with reference feedforward. 
% The frigate's dynamics is modeled using the Norrbin (1963) nonlinear 
% model, providing realistic simulation of naval vessel behavior under 
% various control scenarios.
%
% Dependencies:      
%   frigate.m   - Frigate dynamics
%   refModel.m  - Reference model for autopilot systems
%
% Author:     Thor I. Fossen
% Date:       2024-04-20
% Revisions:

clearvars;

%% USER INPUTS
T_final = 300;	                % Final simulation time (s)
h = 0.05;                       % Sample time (s)

%% AUTOPILOTS PARAMETERS
U = 10;                         % Speed (m/s)
delta_max = deg2rad(30);        % Max rudder angle (deg)

% Reference model parameters
wn_d = 0.1;                     % Natural frequency (rad/s)
zeta_d = 1.0;                   % Relative damping factor (-)
r_max = deg2rad(5.0);           % Maximum turning rate (rad/s)

% PID heading autopilot 
T = 27.0;                       % Nomoto gains at U = 9 m/s
K = 0.18;
wn = 0.5;                       % Closed-loop natural frequency (rad/s)
zeta = 1.0;                     % Closed loop relative damping factor (-)

Kp = (T/K) * wn^2;                    
Kd = (T/K) * (2 * zeta * wn - 1/T);
Td = Kd / Kp; 
Ti = 10 / wn;

% Initial states
psi = 0;                        % Heading angle (rad)
r = 0;                          % Yaw rate (rad/s)
delta = 0;                      % Rudder angle (rad)
psi_int = 0;                    % Integral state
psi_d = psi;                    % Reference model states
r_d = r;
a_d = 0;

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% Display simulation options
displayControlMethod();

%% MAIN LOOP
simdata = zeros(nTimeSteps,6);        % Preallocate table 

for i = 1:nTimeSteps

    % Reference model, step input
    psi_ref = deg2rad(30);
    if t(i) > 100
        psi_ref = deg2rad(70);
    end

   % Reference model propagation
    [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
            zeta_d, wn_d, h, 1);

    delta_c = (T/K) * a_d + (1/K) * r_d -...
        Kp * (ssa( psi-psi_d) +...
        Td * (r - r_d) + (1/Ti) * psi_int );

    % Rudder saturation
    delta_c = sat(delta_c, delta_max);

    % Frigate dynamics
    [psi_dot, r_dot, delta_dot] = frigate(r, U, delta, delta_c);

    % Store data for presentation
    simdata(i,:) = [psi, r, delta, delta_c, psi_d, r_d];

    % Euler's integration method (k+1)
    psi = euler2(psi_dot, psi, h);
    r = euler2(r_dot, r, h);
    delta = euler2(delta_dot, delta, h);
    psi_int = psi_int + h * ssa( psi - psi_d );

end

%% PLOTS
% simdata(i,:) = [psi, r, delta, delta_c, psi_d, r_d] 
psi     = rad2deg(simdata(:,1)); 
r       = rad2deg(simdata(:,2));          
delta   = rad2deg(simdata(:,3));   
delta_c = rad2deg(simdata(:,4));
psi_d   = rad2deg(simdata(:,5));
r_d     = rad2deg(simdata(:,6));

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

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', '100 m', ...
    'Max speed', '12 m/s', ...
    'Max rudder angle', '30 deg'};
displayVehicleData('Frigate', vesselData, 'frigate.jpg', 2);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Frigate');
    disp('Norrbin (1963) nonlinear model');    
    disp('Heading autopilot: PID control law with reference feedforward')
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end
