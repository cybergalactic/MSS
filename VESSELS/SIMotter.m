% SIMotter    User editable script for simulation of the Otter USV under
%             feedback control
%
% Calls:      otter.m
%
% Author:    Thor I. Fossen
% Date:      2021-04-25
% Revisions: 

clearvars;

%% USER INPUTS
h  = 0.02;        % sampling time [s]
N  = 1000;		  % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);	   

% propeller revolutions (rps)
n = [80 60]';          % n = [ n_left n_right ]' 

% Load condition
mp = 25;               % payload mass (kg), max value 45 kg
rp = [0 0 -0.35]';     % location of payload (m)

% Current
V_c = 0;               % current speed (m/s)
beta_c = 30 * pi/180;  % current direction (rad)

%% MAIN LOOP
simdata = zeros(N+1,13);                   % table for simulation data

for i=1:N+1
   t = (i-1) * h;                          % time (s)             

   % store simulation data in a table   
   simdata(i,:) = [t x'];    
   
   % Euler integration (k+1)
   x = x + h * otter(x,n,mp, rp, V_c,beta_c);
end

%% PLOTS
t    = simdata(:,1); 
nu   = simdata(:,2:7); 
eta  = simdata(:,8:13); 

clf
figure(1)
subplot(611),plot(t,nu(:,1),'linewidt',2)
xlabel('time (s)'),title('Surge velocity (m/s)'),grid
subplot(612),plot(t,nu(:,2),'linewidt',2)
xlabel('time (s)'),title('Sway velocity (m/s)'),grid
subplot(613),plot(t,nu(:,3),'linewidt',2)
xlabel('time (s)'),title('Heave velocity (m/s)'),grid
subplot(614),plot(t,(180/pi)*nu(:,4),'linewidt',2)
xlabel('time (s)'),title('Roll rate (deg/s)'),grid
subplot(615),plot(t,(180/pi)*nu(:,5),'linewidt',2)
xlabel('time (s)'),title('Pitch rate (deg/s)'),grid
subplot(616),plot(t,(180/pi)*nu(:,6),'linewidt',2)
xlabel('time (s)'),title('Yaw rate (deg/s)'),grid

figure(2)
subplot(611),plot(eta(:,2),eta(:,1),'linewidt',2); 
title('xy plot (m)'),grid
subplot(612),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2),'linewidt',2);
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(613),plot(t,eta(:,3),'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
subplot(614),plot(t,(180/pi)*eta(:,4),'linewidt',2)
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(615),plot(t,(180/pi)*eta(:,5),'linewidt',2)
xlabel('time (s)'),title('pitch angle (deg)'),grid
subplot(616),plot(t,(180/pi)*eta(:,6),'linewidt',2)
xlabel('time (s)'),title('yaw angle (deg)'),grid