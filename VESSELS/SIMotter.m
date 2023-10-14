% SIMotter    User editable script for simulation of the Otter USV under
%             feedback control
%
% Calls:      otter.m
%
% Author:    Thor I. Fossen
% Date:      2021-04-25
% Revisions: 2023-10-14 Added a heading autopilot and reference model

clearvars;

%% USER INPUTS
h  = 0.02;        % sampling time [s]
N  = 1000;		  % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);	   

% propeller revolutions (rps)
n = [80 60]';             % n = [ n_left n_right ]' 

% Load condition
mp = 25;                  % payload mass (kg), max value 45 kg
rp = [0.05 0 -0.35]';     % location of payload (m)

% Ocean current
V_c = 0.3;                % current speed (m/s)
beta_c = deg2rad(30);     % current direction (rad)

% PID course autopilot (Nomoto gains)
T = 1;
m = 41.4;        % m = T/K
K = T / m;

wn = 1.5;        % pole placement parameters
zeta = 1;

Kp = m * wn^2;
Kd = m * (2*zeta*wn - 1/T);
Td = Kd/Kp; 
Ti = 10/wn;

B = 0.0111 * [1 1                   % input matrix
              0.395 -0.395 ];
Binv = inv(B);

z_psi = 0;                           % integral state

% Reference model
wn_d = 0.5;         % natural frequency
zeta_d = 1.0;       % relative damping factor

psi_d = 0;
r_d = 0;
a_d = 0;

% Propeller dynamics
Tn = 1;                % Propeller time constant      
n = [0 0]';            % n = [ n_left n_right ]'


%% MAIN LOOP
simdata = zeros(N+1,15);                   % table for simulation data

for i=1:N+1

   t = (i-1) * h;                          % time (s)    

   % heading angle setpoints
   if ( t < 20 )
       psi_ref = deg2rad(20);
   else
       chi_ref = deg2rad(0); 
   end

   % Measurements
   r = x(6);
   psi = x(12);

   % Heading autopilot (PID pole placement with reference model)
   tau_X = 100;
   tau_N = (T/K) * a_d + (1/K) * r_d -... 
        Kp * (ssa( psi-psi_d) +...
        Td * (r - r_d) + (1/Ti) * z_psi );
   u = Binv * [tau_X tau_N]';
   n_c = sign(u) .* sqrt( abs(u) );   
   
   % n_c = [80 60]';          % n = [ n_left n_right ]' 

   % Reference model jerk
   j_d = wn_d^3 * ssa(psi_ref - psi_d) - (2*zeta_d+1) * wn_d^2 * r_d...
       - (2*zeta_d+1) * wn_d * a_d;

   % store simulation data in a table   
   simdata(i,:) = [t x' r_d psi_d];    
   
   % Forward Euler (k+1)
   x = x + h * otter(x,n,mp, rp, V_c,beta_c);
   n = n - h/Tn * (n - n_c); 
   z_psi = z_psi + h * ssa( psi-psi_d );
   psi_d = psi_d + h * r_d; 
   r_d = r_d + h * a_d;
   a_d = a_d + h * j_d; 
   
end

%% PLOTS
t    = simdata(:,1); 
nu   = simdata(:,2:7); 
eta  = simdata(:,8:13); 
r_d = simdata(:,14); 
psi_d = simdata(:,15); 

clf
figure(1)
subplot(611),plot(t,nu(:,1))
xlabel('time (s)'),title('Surge velocity (m/s)'),grid
subplot(612),plot(t,nu(:,2))
xlabel('time (s)'),title('Sway velocity (m/s)'),grid
subplot(613),plot(t,nu(:,3))
xlabel('time (s)'),title('Heave velocity (m/s)'),grid
subplot(614),plot(t,rad2deg(nu(:,4)))
xlabel('time (s)'),title('Roll rate (deg/s)'),grid
subplot(615),plot(t,rad2deg(nu(:,5)))
xlabel('time (s)'),title('Pitch rate (deg/s)'),grid
subplot(616),plot(t,rad2deg(nu(:,6)),t,rad2deg(r_d))
xlabel('time (s)'),title('Yaw rate (deg/s)'),grid
legend('r','r_d')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(611),plot(eta(:,2),eta(:,1)); 
title('xy plot (m)'),grid
subplot(612),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2));
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(613),plot(t,eta(:,3),'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
subplot(614),plot(t,rad2deg(eta(:,4)))
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(615),plot(t,rad2deg(eta(:,5)))
xlabel('time (s)'),title('pitch angle (deg)'),grid
subplot(616),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d))
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('\psi','\psi_d')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
