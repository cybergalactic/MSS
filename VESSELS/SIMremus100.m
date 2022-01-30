% SIMremus100 User editable script for simulation of the Remus 100 AUV under
%             feedback control
%
% Calls:      remus100.m
%
% Author:     Thor I. Fossen
% Date:       2021-06-28
% Revisions:  2021-12-16 retuning of the autpilots, added integral action

clearvars;

%% USER INPUTS
h  = 0.05;               % sample time (s)
N  = 20000;              % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);
z_int = 0;

% initial setpoints
n = 0;                  % propeller revolution (rpm)
z_d = 0;                % depth (m)
psi_d = 0;              % heading angle (rad)

% ocean current velcoities expressed in NED
Vc = 0.5;                     % speed (m/s)
betaVc = -10 * pi/180;        % direction (rad)

% controller gains
Kp_z = 0.1;                   % depth controller
T_z = 1000;
Kp_theta = 1.0;             
Kd_theta = 1.0;

Kp_psi = 0.5;                 % heading autopilot
Kd_psi = 1; 
   
%% MAIN LOOP
simdata = zeros(N+1,16);                   % table for simulation data
for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % depth controller (succesive-loop closure)    
   if t > 100
       z_ref = 20;
   else
       z_ref = 0;
   end
   z_d = exp(-h/10) * z_d + (1 - exp(-h/10)) * z_ref;  % LP filter
   
   theta_d = Kp_z * ( (x(9) - z_d) + (1/T_z) * z_int );                % PI 
   delta_s = -Kp_theta * ( ssa( x(11) - theta_d ) - Kd_theta * x(5) ); % PD
   
   % heading controller
   if t > 150
      psi_ref = -40 * pi/180;
   else
       psi_ref = 0;
   end   
   psi_d = exp(-h/5) * psi_d + (1 - exp(-h/5)) * psi_ref;  % LP filter
   
   delta_r = -Kp_psi * ssa( x(12) - psi_d ) - Kd_psi * x(6);           % PD

   % propeller revolution (rpm)
   if (n < 1500)
       n = n + 1;
   end
   
   % control inputs
   ui = [delta_r delta_s n]';
   
   % store simulation data in a table 
   simdata(i,:) = [t x' ui'];   
   
   % Euler integration (k+1)
   z_int = z_int + h * (x(9) - z_d);
   x = x + h * remus100(x,ui,Vc,betaVc);    
   
end

%% PLOTS
t   = simdata(:,1);  
nu  = simdata(:,2:7);  
eta = simdata(:,8:13);  
u   = simdata(:,14:16); 

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
subplot(612),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2+nu(:,3).^2),'linewidt',2);
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(613),plot(t,eta(:,3),'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
subplot(614),plot(t,(180/pi)*eta(:,4),'linewidt',2)
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(615),plot(t,(180/pi)*eta(:,5),'linewidt',2)
xlabel('time (s)'),title('pitch angle (deg)'),grid
subplot(616),plot(t,(180/pi)*eta(:,6),'linewidt',2)
xlabel('time (s)'),title('yaw angle (deg)'),grid

figure(3)
subplot(311),plot(t,(180/pi)*u(:,1),'linewidt',2)
xlabel('time (s)'),title('Rudder \delta_r (deg)'),grid
subplot(312),plot(t,(180/pi)*u(:,2),'linewidt',2)
xlabel('time (s)'),title('Stern planes \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3),'linewidt',2)
xlabel('time (s)'),title('Propeller revolutions n (rpm)'),grid

