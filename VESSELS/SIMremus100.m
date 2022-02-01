% SIMremus100 User editable script for simulation of the Remus 100 AUV under
%             feedback control
%
% Calls:      remus100.m
%
% Author:     Thor I. Fossen
% Date:       2021-06-28
% Revisions:  2022-02-01 redesign of the autopilots

clearvars;

%% USER INPUTS
h  = 0.05;               % sample time (s)
N  = 40000;              % number of samples

% initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1);
z_int = 0;
theta_int = 0;
psi_int = 0;

% setpoints
n = 0;                  % initial propeller revolution (rpm)
n_d = 1500;             % desired propeller revolution, max 1500 rpm
z_d = 0;                % initial depth (m)
z_step = 20;            % step change in depth, max 100 m
psi_d = 0;              % initial heading angle (rad)
psi_step = -20*pi/180;  % step change in heading angle (rad)

% ocean current velcoities expressed in NED
Vc = 0.5;                     % speed (m/s)
betaVc = -10 * pi/180;        % direction (rad)

% depth controller
wn_d_z = 1/500;             % desired natural frequency, reference model
Kp_z = 0.05;                 
T_z = 1000;

wn_b_theta = 5;             % bandwidth, pole placement algorithm 
m55 = 8.6;                  % moment of inertia, pitch
d55 = 11.9;                 % linear damping, pitch
Kp_theta = m55 * wn_b_theta^2;             
Kd_theta = m55 * wn_b_theta - d55;
Ki_theta = Kp_theta * (wn_b_theta/10);

% heading autopilot 
wn_d_psi = 1/20;            % desired natural frequency, reference model 
wn_b_psi = 1;               % bandwidth, pole placement algorithm 
m66 = 8.6;                  % moment of inertia, yaw
d66 = 0.6;                  % linear damping, yaw
Kp_psi = m66 * wn_b_psi^2;                 
Kd_psi = m66 * 2*wn_b_psi - d66; 
Ki_psi = Kp_psi * (wn_b_psi/10);
   
%% MAIN LOOP
simdata = zeros(N+1,19);                   % table for simulation data
for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % depth controller (succesive-loop closure)  
   if (z_step > 100 || z_step < 0)
       error('desired depth must be between 0-100 m')
   end  
   
   if t > 100
       z_ref = z_step;
   else
       z_ref = 0;
   end
   z_d = exp(-h*wn_d_z) * z_d + (1 - exp(-h*wn_d_z)) * z_ref;  % LP filter
   
   theta_d = Kp_z * ( (x(9) - z_d) + (1/T_z) * z_int );                % PI 
   delta_s = -Kp_theta * ssa( x(11) - theta_d ) - Kd_theta * x(5)...
       - Ki_theta * theta_int;                                        % PID
   
   if abs(delta_s) > 70 * (pi/180)
       delta_s = sign(delta_s) * 70 * (pi/180);
   end
   
   % heading controller
   if t > 150
      psi_ref = psi_step;
   else
       psi_ref = 0;
   end   
   
   % LP filter
   psi_d = exp(-h*wn_d_psi) * psi_d + (1 - exp(-h*wn_d_psi)) * psi_ref;     
   
   delta_r = -Kp_psi * ssa( x(12) - psi_d ) - Kd_psi * x(6)...
        -Ki_psi * psi_int;                                           % PID 

   % propeller revolution (rpm)
   if (n < n_d)
       n = n + 1;
   end
   
   % control inputs
   ui = [delta_r delta_s n]';
   
   % store simulation data in a table 
   simdata(i,:) = [t x' ui' z_d theta_d psi_d ];   
   
   % Euler integration (k+1)
   z_int = z_int + h * (x(9) - z_d);
   theta_int = theta_int + h * ssa( x(11) - theta_d );
   psi_int = psi_int + h * ssa( x(12) - psi_d );   
   x = x + h * remus100(x,ui,Vc,betaVc);    
   
end

%% PLOTS
t       = simdata(:,1);  
nu      = simdata(:,2:7);  
eta     = simdata(:,8:13);  
u       = simdata(:,14:16); 
z_d     = simdata(:,17); 
theta_d = simdata(:,18); 
psi_d   = simdata(:,19); 

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
subplot(613),plot(t,eta(:,3),t,z_d,'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
legend('true','desired')
subplot(614),plot(t,(180/pi)*eta(:,4),'linewidt',2)
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(615),plot(t,(180/pi)*eta(:,5),t,(180/pi)*theta_d,'linewidt',2)
xlabel('time (s)'),title('pitch angle (deg)'),grid
legend('true','desired')
subplot(616),plot(t,(180/pi)*eta(:,6),t,(180/pi)*psi_d,'linewidt',2)
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('true','desired')

figure(3)
subplot(311),plot(t,(180/pi)*u(:,1),'linewidt',2)
xlabel('time (s)'),title('Rudder \delta_r (deg)'),grid
subplot(312),plot(t,(180/pi)*u(:,2),'linewidt',2)
xlabel('time (s)'),title('Stern planes \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3),'linewidt',2)
xlabel('time (s)'),title('Propeller revolutions n (rpm)'),grid

