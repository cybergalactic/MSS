function SIMmariner()
% SIMmariner is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script simulates the dynamic behavior of a mariner class vessel 
% under PID heading control. This simulation utilizes a PID control 
% strategy to maintain heading and demonstrates the vessel’s  response 
% visually through an animation of its trajectory in the North-East 
% coordinate plane. An extension to waypoint path-following control using 
% a course autopilot fpr turning is included in Simulink demo library.
%
% Dependencies:
%   mariner.m - Vessel dynamics.  
%
% Simulink Models:
%   demoMarinerPathFollowingCourseControl.slx
%
% Author:     Thor I. Fossen
% Date:       2018-07-21
% Revisions:
%   2024-03-27 : Added animation of the ship's North-East positions.
%   2024-04-19 : Enhanced compatibility with GNU Octave.

close all;
clear animateShip  % clear the persistent animation variables

t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)

Kp = 1;      % controller proportional gain
Td = 10;     % controller derivative time (s)
Ti = 100;    % controller integral time (s)

% Initial states
x = zeros(7,1);     % x = [ u v r x y psi delta ]'
z_psi = 0;          % integral state

%% MAIN LOOP
N = round(t_f/h);                    % number of samples
simdata = zeros(N+1,length(x)+2);    % memory allocation

for i=1:N+1

    time = (i-1) * h;                % simulation time in seconds

    % Measurements
    r   = x(3) + 0.001 * randn;
    psi = x(6) + 0.01 * randn;
    
    % PID control system
    psi_ref = deg2rad(5);            % desired heading
    delta = -Kp * ( ssa(psi - psi_ref) + Td * r + 1/Ti * z_psi );  

    % Ship dynamics
    [xdot,U] = mariner(x,delta);     
    
    % Store data for presentation
    simdata(i,:) = [time,x',U]; 
    
    % Numerical integration
    x = x + h * xdot;                             % Euler's method
    z_psi = z_psi + h * ssa(psi  - psi_ref);

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

% Simdata(i,:) = [t, x', U]
t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));
delta = rad2deg(simdata(:,8));
U     = simdata(:,9);

% Plot and animation of the North-East positions
figure(1)
shipSize = 1.0;
set(gcf, 'Position', [1, 1, 0.4*scrSz(3), scrSz(4)]);
animateShip(x,y,shipSize,'b-',1);

figure(2)
subplot(221)
plot(t,r)
xlabel('Time (s)')
title('yaw rate r (deg/s)')
grid
subplot(222)
plot(t,U)
xlabel('Time (s)')
title('speed U (m/s)')
grid
subplot(223)
plot(t,psi,[0,t(end)],rad2deg([psi_ref psi_ref]))
xlabel('Time (s)')
title('yaw angle \psi (deg)')
grid
subplot(224)
plot(t,delta)
xlabel('Time (s)')
title('rudder angle \delta (deg)')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

end

