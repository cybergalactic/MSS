function SIMdsrv()
% SIMdsrv is compatibel with MATLAB and GNU Octave (www.octave.org). This
% script simulates the Naval Postraduate School's Deep Submergence Rescue 
% vehicle (DSRV) under depth-changing maneuvers. The depth autopilot is 
% designed using using succesive-loop closure and PID methods.
%
% Dependencies:
%   DSRV.m    - DSRV dynamics
%
% Author:     Thor I. Fossen
% Date:       2024-04-20
% Revisions:
%   2024-04-19 : Enhanced compatibility with GNU Octave.
%   2024-07-10: Improved numerical accuracy by replacing Euler's method with RK4

clearvars;

%% USER INPUTS
T_final = 1000;	            % Final simulation time (s)
h  = 0.05;                  % Sample time (s)

%% DEPTH AUTOPILOTS PARAMETERS
z_d = 0;                    % Initial depth (m), reference model
z_step = 200;               % Step change in depth

% Autopilot integral states
z_int = 0;                  % Depth
theta_int = 0;              % Pitch angle

% Stern rudder
T_delta = 0.1;              % Time constant (s)
delta_max = deg2rad(30);    % Max stern plane angle (deg)

% Depth controller (suceessive-loop closure)
wn_d_z = 0.01;              % Desired natural frequency, reference model
wnz = 1 * wn_d_z;           % Desired natural frequency (heave)
Kp_z = 0.01;                % Proportional gain (heave)
T_z = 100;                  % Integral time constant (heave)
Kp_theta = 1.0;             % Proportional gain (pitch)
Ki_theta = 0.1;             % Integral gain (pitch)

% Initial states
x = zeros(5,1);     % x = [ w q x z theta ]'
delta_s = 0;

% Display simulation options
displayControlMethod();

%% MAIN LOOP
t = 0:h:T_final;                     % Time vector
simdata = zeros(length(t),8);        % Preallocate table 

for i = 1:length(t)

    % Measurements
    w     = x(1) + 0.001 * randn;
    q     = x(2) + 0.001 * randn;
    z     = x(4) + 0.001 * randn;
    theta = x(5) + 0.001 * randn;

    % Depth command, z_ref
    if t(i) > 300
        z_ref = z_step;
    else
        z_ref = 10;
    end

    % LP filtering the depth command
    z_d = exp(-h*wnz) * z_d + (1 - exp(-h*wnz)) * z_ref;

    % Depth autopilot using the stern planes (succesive-loop closure)
    theta_d = Kp_z * ( (z - z_d) + (1/T_z) * z_int );    
    delta_PID = -Kp_theta * ssa( theta - theta_d ) - Ki_theta * theta_int;
    delta_c = -delta_PID;
    delta_c = sat(delta_c, delta_max);   % Amplitude saturation

    % Store data for presentation
    simdata(i,:) = [x', delta_s, theta_d, z_d];

    % RK4 method (k+1)
    x = rk4(@DSRV, h, x, delta_s);  % DSRV dynamics

    % Euler's integration method (k+1)
    delta_s = delta_s + h * (delta_c - delta_s) / T_delta;
    z_int = z_int + h * ( z - z_d );
    theta_int = theta_int + h * ssa( theta - theta_d );

end

%% PLOTS
% simdata(i,:) = [ w q x z theta delta_s theta_d z_d ]   
w       = simdata(:,1); 
q       = rad2deg(simdata(:,2));          
x       = simdata(:,3);   
z       = simdata(:,4);
theta   = rad2deg(simdata(:,5));
delta_s = rad2deg(simdata(:,6));
theta_d = rad2deg(simdata(:,7));
z_d     = simdata(:,8);

figure(1); clf;
plot(t,delta_s)
xlabel('Time (s)')
title('Stern rudder angle \delta_s (deg)')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

figure(2); clf;
subplot(221)
plot(t,w)
xlabel('Time (s)')
title('Heave velocity (m/s)')
grid
subplot(222)
plot(t,q)
xlabel('Time (s)')
title('Pitch rate q (deg/s)')
grid
subplot(223)
plot(t,theta,t,theta_d)
xlabel('Time (s)')
title('Pitch angle \theta (deg)')
grid
legend('True','Desired')
subplot(224)
plot(t,z,t,z_d)
xlabel('Time (s)')
title('Depth z (m)')
legend('True','Desired')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', '5.0 m',...  
    'Mass', '2300 kg',...
    'Maximum speed', '4.11 m/s',...
    'Max rudder angle', '30 deg'};
displayVehicleData('Deep Submergence Rescue Vehicle (DSRV)', ...
    vesselData, 'DSRV.jpg', 3);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Deep Submergence Rescue Vehicle (DSRV)');
    disp('Depth autopilot: Succesive-loop closure')
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end

