% exOtter is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates the Maritime Robotics Otter Uncrewed Surface 
% Vehicle (USV), which is controlled by two propellers. SOG, COG, and 
% course rate are estimated using an EKF, and the course autopilot is 
% utilizing a PID pole-placement algorithm.
%
% Author:    Thor I. Fossen
% Date:      2019-07-18
% Revisions: 
%   2021-05-25 Added EKF for COG/SOG/course rate and course autopilot
%   2024-07-10: Updated to match SIMotter.m and use RK4

clearvars;
clear EKF_5states   % Clear persistent states in EKF_5states.m

%% Simulation data
fHz = 50;
h  = 1/fHz;         % Sampling time [s]
Z = 10;             % GNSS measurement frequency (10 times slower)
T_final = 40;	    % Final simulation time [s]

% Initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1); x(1) = 1;	   

% EKF covariance matrices and initial states
Qd = 500 * diag([ 1000 1000 ]);
Rd = 0.00000001 * diag([ 1 1 ]);
x_hat = zeros(5,1);

% Otter USV input matrix
[~,~,M, B_prop] = otter();
Binv = inv(B_prop);              % Invert input matrix for control allocation

% PID heading autopilot parameters (Nomoto model: M(6,6) = T/K)
T = 1;                           % Nomoto time constant
K = T / M(6,6);                  % Nomoto gain constant

wn = 1.5;                        % Closed-loop natural frequency (rad/s)
zeta = 1.0;                      % Closed-loop relative damping factor (-)

Kp = M(6,6) * wn^2;                     % Proportional gain
Kd = M(6,6) * (2 * zeta * wn - 1/T);    % Derivative gain
Td = Kd / Kp;                           % Derivative time constant
Ti = 10 / wn;                           % Integral time constant

z = 0;                              % Integral state

% Reference model parameters
wn_d = 1.0;                      % Natural frequency (rad/s)
zeta_d = 1.0;                    % Relative damping factor (-)
r_max = deg2rad(10.0);           % Maximum turning rate (rad/s)

chi_d = 0;
omega_d = 0;
a_d = 0;

% Propeller dynamics
Tn = 0.1;              % Propeller time constant      
n = [0 0]';            % n = [ n_left n_right ]'

% Load condition
mp = 25;                         % Payload mass (kg), maximum value 45 kg
rp = [0.05 0 -0.35]';            % Location of payload (m)

% Ocean current
V_c = 0;               % Current speed (m/s)
beta_c = 30 * pi/180;  % Current direction (rad)

%% MAIN LOOP
t = 0:h:T_final;                 % Time vector
simdata = zeros(length(t),19);   % Table for simulation data
for i=1:length(t)

   % Course angle setpoints
   if ( t(i) < 20 )
       chi_ref = deg2rad(20);
   else
       chi_ref = deg2rad(0); 
   end
   
   % Course autopilot (PID pole placement with reference model)
   tau_X = 100;
   tau_N = (T/K) * a_d + (1/K) * omega_d -... 
        Kp * ( ssa( x_hat(4)-chi_d ) +...
        Td * (x_hat(5) - omega_d) + (1/Ti) * z );
   u = Binv * [tau_X tau_N]';
   n_c = sign(u) .* sqrt( abs(u) );   
   
   % n_c = [80 60]';          % n = [ n_left n_right ]' 
   
   % Reference model propagation
   [chi_d, omega_d, a_d] = refModel(chi_d, omega_d, a_d, chi_ref, ...
      r_max, zeta_d, wn_d, h, 1);
  
   % Store simulation data in a table   
   simdata(i,:) = [x' x_hat' n'];    
   
   % 5-state EKF:  x_hat[k+1] = [x_N, v_N, U, chi, omega]'
   x_N = x(7); y_E = x(8);      % GNSS measurements
   x_hat = EKF_5states(x_N,y_E,h,Z,'NED',Qd,Rd);  
         
   % RK4 method x(k+1)
   x = rk4(@otter, h, x, n,mp,rp,V_c,beta_c);
   
   % Euler's method (k+1)
   z = z + h * ssa( x_hat(4)-chi_d );
   n = n - h / Tn * (n - n_c);  

end

%% PLOTS
nu   = simdata(:,1:6); 
eta  = simdata(:,7:12); 

u = nu(:,1);
v = nu(:,2);
U = sqrt(u.^2 + v.^2);           % Speed
psi = eta(:,6);                  % Heading angle
beta_c = ssa( atan2(v,u) );      % Crab angle
chi = psi + beta_c;              % Course angle

x_hat = simdata(:,13);
y_hat = simdata(:,14);
U_hat = simdata(:,15);
chi_hat = simdata(:,16);
omega_hat = simdata(:,17);

n1 = simdata(:,18);
n2 = simdata(:,19);

clf
figure(1)
subplot(611),plot(t,nu(:,1))
xlabel('time (s)'),title('Surge velocity (m/s)'),legend('True'),grid
subplot(612),plot(t,nu(:,2))
xlabel('time (s)'),title('Sway velocity (m/s)'),legend('True'),grid
subplot(613),plot(t,nu(:,3))
xlabel('time (s)'),title('Heave velocity (m/s)'),legend('True'),grid
subplot(614),plot(t,(180/pi)*nu(:,4))
xlabel('time (s)'),title('Roll rate (deg/s)'),legend('True'),grid
subplot(615),plot(t,(180/pi)*nu(:,5))
xlabel('time (s)'),title('Pitch rate (deg/s)'),legend('True'),grid
subplot(616),plot(t,(180/pi)*nu(:,6))
xlabel('time (s)'),title('Yaw rate (deg/s)'),legend('True'),grid
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

figure(2)
subplot(611),plot(eta(:,2),eta(:,1),y_hat,x_hat); 
title('xy plot (m)'),legend('True','Estimate'),grid
subplot(612),plot(t,U,t,U_hat);
xlabel('time (s)'),title('speed (m/s)'),legend('True','Estimate'),grid
subplot(613),plot(t,eta(:,3))
xlabel('time (s)'),title('heave z position (m)'),legend('True'),grid
subplot(614),plot(t,(180/pi)*eta(:,4))
xlabel('time (s)'),title('roll angle (deg)'),legend('True'),grid
subplot(615),plot(t,(180/pi)*eta(:,5))
xlabel('time (s)'),title('pitch angle (deg)'),legend('True'),grid
subplot(616),plot(t,(180/pi)*psi)
xlabel('time (s)'),title('yaw angle (deg)'),
legend('True'),grid
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

figure(3)
subplot(311),plot(t,(180/pi)*chi,t,(180/pi)*chi_hat)
xlabel('time (s)'),title('course angle (deg)')
legend('True','Estimate'),grid
subplot(312),plot(t,(180/pi)*omega_hat,'r')
xlabel('time (s)'),title('course rate (deg(/s)')
legend('Estimate'),grid
subplot(313),plot(t,n1,'g',t,n2,'k'),grid
xlabel('time (s)')
legend('n_1 left propeller','n_2 right propeller')
title('Propeller revolutions (rad/s)')
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
