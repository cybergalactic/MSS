echo on
% SIMmariner  User editable script for simulation of the mariner class 
% vessel under PID heading control.
%
% Calls:      mariner.m
%             euler2.m
%
% Author:     Thor I. Fossen
% Date:       2018-07-21
% Revisions:  2024-03-27 Added animation of the ship North-East positions

clear animateShip  % clear the persistent animation variables
clearvars;

echo off 
disp('Simulating mariner.m under PD-control with psi_ref = 5 (deg)...')

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
    x = euler2(xdot,x,h);                           % Euler's method
    z_psi = z_psi + h * ssa(psi  - psi_ref);

end

%% PLOTS
t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));
delta = rad2deg(simdata(:,8));
U     = simdata(:,9);

figNo = 1;
shipSize = 1.0;
animateShip(x,y,shipSize,'b-',figNo); % animation of the North-East positions

figure(2)
subplot(221)
plot(t,r)
xlabel('time (s)')
title('yaw rate r (deg/s)')
grid
subplot(222)
plot(t,U)
xlabel('time (s)')
title('speed U (m/s)')
grid
subplot(223)
plot(t,psi,[0,t(end)],rad2deg([psi_ref psi_ref]))
xlabel('time (s)')
title('yaw angle \psi (deg)')
grid
subplot(224)
plot(t,delta)
xlabel('time (s)')
title('rudder angle \delta (deg)')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

