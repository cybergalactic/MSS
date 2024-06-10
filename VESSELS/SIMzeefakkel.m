function SIMzeefakkel()
% SIMzeefakkel is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates a heading-controlled ship, the Zeefakkel with a 
% length of 45 meters, employing a PID control strategy with reference 
% feedforward. The ship dynamics is modeled using the Norrbin (1963) 
% nonlinear model, providing realistic simulation of the ship turning 
% behavior under various control scenarios.
%
% Dependencies:      
%   zeefakkel.m - Craft dynamics
%   refModel.m  - Reference model for autopilot systems
%
% Author:     Thor I. Fossen
% Date:       2024-06-10
% Revisions:
%   None

close all;
clearvars;

%% USER INPUTS
h  = 0.05;                      % Sample time (s)
N  = 10000;                     % Number of samples

%% AUTOPILOTS PARAMETERS
U = 6;                          % Speed (m/s)
delta_max = deg2rad(30);        % Max rudder angle (deg)

% Reference model parameters
wn_d = 0.1;                     % Natural frequency (rad/s)
zeta_d = 1.0;                   % Relative damping factor (-)
r_max = deg2rad(5.0);           % Maximum turning rate (rad/s)

% PID heading autopilot 
T = 31.0;                       % Nomoto gains at U = 5 m/s
K = 0.5;
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

% Display simulation options
displayControlMethod();

%% MAIN LOOP
simdata = zeros(N+1,7);              % Memory allocation

for i = 1:N+1

    t= (i-1) * h;                    % Simulation time in seconds

    % Reference model, step input
    psi_ref = deg2rad(-100);
    if t > 100
        psi_ref = deg2rad(60);
    end

   % Reference model propagation
    [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
            zeta_d, wn_d, h, 1);

    delta_c = (T/K) * a_d + (1/K) * r_d -...
        Kp * (ssa( psi-psi_d) +...
        Td * (r - r_d) + (1/Ti) * psi_int );

    delta_c = sat(delta_c, delta_max);   % Amplitude saturation

    % Craft dynamics
    [psi_dot, r_dot, delta_dot] = zeefakkel(r, U, delta, delta_c);

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
r_d     = rad2deg(simdata(:,7));

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
    'Length', '45 m', ...
    'Beam', '8 m', ...
    'Max speed', '7 m/s', ...
    'Max rudder angle', '30 deg'};
displayVehicleData('Zeefakkel', vesselData, 'zeefakkel.jpg', 2);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Zeefakkel');
    disp('Norrbin (1963) nonlinear model');    
    disp('Heading autopilot: PID control law with reference feedforward')
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end
