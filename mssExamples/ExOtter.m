% ExOtter 
%
% Author:    Thor I. Fossen
% Date:      2019-07-18
% Revisions: 

clearvars;

%% USER INPUTS
h  = 0.02;        % sampling time [s]
N  = 1000;		  % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);	   

% propeller revolutions (rad/s(
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
subplot(611),plot(t,nu(:,1))
xlabel('time (s)'),title('Surge velocity'),grid
subplot(612),plot(t,nu(:,2))
xlabel('time (s)'),title('Sway velocity'),grid
subplot(613),plot(t,nu(:,3))
xlabel('time (s)'),title('Heave velocity'),grid
subplot(614),plot(t,(180/pi)*nu(:,4))
xlabel('time (s)'),title('Roll rate'),grid
subplot(615),plot(t,(180/pi)*nu(:,5))
xlabel('time (s)'),title('Pitch rate'),grid
subplot(616),plot(t,(180/pi)*nu(:,6))
xlabel('time (s)'),title('Yaw rate'),grid

figure(2)
subplot(611),plot(eta(:,2),eta(:,1)); 
title('xy plot'),grid
subplot(612),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2));
xlabel('time (s)'),title('speed'),grid
subplot(613),plot(t,eta(:,3))
xlabel('time (s)'),title('heave z position'),grid
subplot(614),plot(t,(180/pi)*eta(:,4))
xlabel('time (s)'),title('roll angle'),grid
subplot(615),plot(t,(180/pi)*eta(:,5))
xlabel('time (s)'),title('pitch angle'),grid
subplot(616),plot(t,(180/pi)*eta(:,6))
xlabel('time (s)'),title('yaw angle'),grid
