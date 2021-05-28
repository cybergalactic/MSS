% SIMremus100 User editable script for simulation of the Remus 100 AUV under
%             feedback control
%
% Calls:      remus100.m
%
% Author:    Thor I. Fossen
% Date:      2021-06-28
% Revisions: 

clearvars;

%% USER INPUTS
h  = 0.05;               % sample time (s)
N  = 10000;              % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);

% initial setpoints
n = 0;                  % propeller revolution (rpm)
z_d = 0;                % depth (m)
psi_d = 0;              % heading angle (rad)

% ocean current velcoities expressed in BODY
v_current = [0.5 -0.1 0]';

% controller gains
Kp_z = 0.1;                   % depth controller
Kp_theta = 0.2;             
Td_theta = 1.0;

Kp_psi = 0.8;                 % heading autopilot
Kd_psi = 1; 
   
%% MAIN LOOP
simdata = zeros(N+1,16);                   % table for simulation data
for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % Depth controller (succesive loop closure)    
   if t > 100
       z_ref = 50;
   else
       z_ref = 0;
   end
   z_d = exp(-h/10) * z_d + (1 - exp(-h/10)) * z_ref;  % LP filter
   
   theta_d = Kp_z * (x(9) - z_d);   
   delta_s = -Kp_theta * ( ssa( x(11) - theta_d ) - Td_theta * x(5) );
   
   % Heading controller
   if t > 150
      psi_ref = -40 * pi/180;
   else
       psi_ref = 0;
   end   
   psi_d = exp(-h/5) * psi_d + (1 - exp(-h/5)) * psi_d;  % LP filter
   
   delta_r = -Kp_psi * ssa( x(12) - psi_ref ) - Kd_psi * x(6);

   % propeller revolution (rpm)
   if (n < 1500)
       n = n + 1;
   end
   
   % control inputs
   ui = [delta_r delta_s n]';
   
   % store simulation data in a table 
   simdata(i,:) = [t x' ui'];   
   
   % Euler integration (k+1)
   x = x + h * remus100(x,ui,v_current);                	   
   
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

