% SIMremus100 User editable script for simulation of the Remus 100 AUV
%             (remus100.m) under feedback control (depth and heading control 
%             when exposed to ocean currents). Both the Euler angle and unit
%             quaternion representations of the Remus 100 model can be used.
%
% Calls:      remus100.m
%
% Simulink:   demoAUVdepthHeadingControl.slx
%
% Author:     Thor I. Fossen
% Date:       2021-06-28
% Revisions:  2022-02-01 redesign of the autopilots
%             2022-05-06 retuning for new version of remus100.m
%             2022-05-08 added compability for unit quaternions

clearvars;

%% USER INPUTS
h  = 0.05;               % sample time (s)
N  = 10000;              % number of samples
param = 1;               % 0 = Euler angles, 1 = unit quaternions

% initial state vector
phi = 0; theta = 0; psi = 0;

if (param == 0)         % x = [ u v w p q r x y z phi theta psi ]'   
    x = [zeros(9,1); phi; theta; psi];
else                    % x = [ u v w p q r x y z eta eps1 eps2 eps3 ]'
    quat = euler2q(phi,theta,psi);
    x = [zeros(9,1); quat];
end

% integral states (autopilots)
z_int = 0;
theta_int = 0;
psi_int = 0;

% setpoints (autopilots)
n = 0;                  % initial propeller revolution (rpm)
n_d = 1525;             % desired propeller revolution, max 1525 rpm
z_d = 0;                % initial depth (m)
z_step = 30;            % step change in depth, max 100 m
psi_d = 0;              % initial heading angle (rad)
psi_step = -60*pi/180;  % step change in heading angle (rad)

% ocean current velcoities expressed in NED
Vc = 0.5;                   % speed (m/s)
betaVc = 170 * pi/180;      % direction (rad)

% depth controller
wn_d_z = 1/20;              % desired natural frequency, reference model
Kp_z = 0.1;                 % proportional gain
T_z = 100;                  % integral gain
Kp_theta = 2;            
Kd_theta = 3;
Ki_theta = 0.1;

% heading autopilot 
wn_d_psi = 1/5;             % desired natural frequency, reference model 
wn_b_psi = 1;               % bandwidth, pole placement algorithm 
m66 = 7.5;                  % moment of inertia, yaw
Kp_psi = m66 * wn_b_psi^2;                 
Kd_psi = m66 * 2*wn_b_psi; 
Ki_psi = Kp_psi * (wn_b_psi/10);
   
%% MAIN LOOP
if (param == 0)
    disp('...simulating the Remus 100 AUV using Euler angles (12 states)')
else
    disp('...simulating the Remus 100 AUV using unit quaternions (13 states)') 
end

simdata = zeros(N+1,length(x)+7); % allocate empty table for simulation data

for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % measurements
   q = x(5);                % pitch rate
   r = x(6);                % yaw rate
   z = x(9);                % z-position (depth)
   
   if (param==0)
         phi = x(10); theta = x(11); psi = x(12);   % Euler angles
   else
         [phi,theta,psi] = q2euler(x(10:13));  % quaternion to Euler angles
   end
   
   % depth controller (succesive-loop closure)  
   if (z_step > 100 || z_step < 0)
       error('desired depth must be between 0-100 m')
   end  
   
   if t > 200
       z_ref = z_step;
   else
       z_ref = 0;
   end
   z_d = exp(-h*wn_d_z) * z_d + (1 - exp(-h*wn_d_z)) * z_ref;  % LP filter
   
   theta_d = Kp_z * ( (z - z_d) + (1/T_z) * z_int );               % PI 
   delta_s = -Kp_theta * ssa( theta - theta_d ) - Kd_theta * q...
       - Ki_theta * theta_int;                                     % PID
   
   % PID heading controller
   if t > 100
      psi_ref = psi_step;
   else
       psi_ref = 0;
   end   
   
   % LP filter
   psi_d = exp(-h*wn_d_psi) * psi_d + (1 - exp(-h*wn_d_psi)) * psi_ref;     
   
   delta_r = -Kp_psi * ssa( psi - psi_d ) - Kd_psi * r - Ki_psi * psi_int;                                           % PID 

   % propeller revolution (rpm)
   if (n < n_d)
       n = n + 1;
   end
   
   % control inputs 
   max_ui = [30*pi/180 30*pi/180  1525]';   % rad, rad, rpm
   if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
   if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
   if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end
    
   ui = [delta_r delta_s n]';
   
   % store simulation data in a table 
   simdata(i,:) = [t z_d theta_d psi_d ui' x'];   
   
   % Propagate the vehicle dynamics (k+1)
   xdot = remus100(x,ui,Vc,betaVc);
   
   if (param == 0)     
       x = x + h * xdot;                       % Euler's integration method
   else
       x(1:9) = x(1:9) + h * xdot(1:9);        % Euler's integration method
       
       quat = x(10:13);                        % unit quaternion
       quat = expm(Tquat(x(4:6)) * h) * quat;  % exact discretization
       x(10:13) = quat/norm(quat);             % normalization
   end
   
   % Euler's integration method (k+1)
   z_int = z_int + h * ( z - z_d );
   theta_int = theta_int + h * ssa( theta - theta_d );
   psi_int = psi_int + h * ssa( psi - psi_d );   
   
end

%% PLOTS
t       = simdata(:,1);         % simdata = [t z_d theta_d psi_d ui' x']
z_d     = simdata(:,2); 
theta_d = simdata(:,3); 
psi_d   = simdata(:,4); 
u       = simdata(:,5:7); 
nu      = simdata(:,8:13);

if (param==0)                   % Euler angle representation
    eta = simdata(:,14:19);
else                            % Transform unit quaternions to Euler angles
    quaternion = simdata(:,17:20);
    for i = 1:N+1
        [phi(i,1),theta(i,1),psi(i,1)] = q2euler(quaternion(i,:));
    end
    eta = [simdata(:,14:16) phi theta psi];
    
    figure(4)
    plot(t,quaternion,'linewidth',2);
    xlabel('time (s)'),title('Unit quaternion'),grid
    legend('\eta','\epsilon_1','\epsilon_2','\epsilon_3')
end

figure(1); 
subplot(611),plot(t,nu(:,1))
xlabel('time (s)'),title('Surge velocity (m/s)'),grid
subplot(612),plot(t,nu(:,2))
xlabel('time (s)'),title('Sway velocity (m/s)'),grid
subplot(613),plot(t,nu(:,3))
xlabel('time (s)'),title('Heave velocity (m/s)'),grid
subplot(614),plot(t,(180/pi)*nu(:,4))
xlabel('time (s)'),title('Roll rate (deg/s)'),grid
subplot(615),plot(t,(180/pi)*nu(:,5))
xlabel('time (s)'),title('Pitch rate (deg/s)'),grid
subplot(616),plot(t,(180/pi)*nu(:,6))
xlabel('time (s)'),title('Yaw rate (deg/s)'),grid

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

figure(2); 
subplot(611),plot(eta(:,2),eta(:,1)); 
title('xy plot (m)'),grid
subplot(612),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2+nu(:,3).^2));
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(613),plot(t,eta(:,3),t,z_d)
xlabel('time (s)'),title('heave position (m)'),grid
legend('true','desired')
subplot(614),plot(t,(180/pi)*eta(:,4))
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(615),plot(t,(180/pi)*eta(:,5),t,(180/pi)*theta_d)
xlabel('time (s)'),title('pitch angle (deg)'),grid
legend('true','desired')
subplot(616),plot(t,(180/pi)*eta(:,6),t,(180/pi)*psi_d)
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('true','desired')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

figure(3); 
subplot(311),plot(t,(180/pi)*u(:,1))
xlabel('time (s)'),title('Rudder \delta_r (deg)'),grid
subplot(312),plot(t,(180/pi)*u(:,2))
xlabel('time (s)'),title('Stern planes \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3))
xlabel('time (s)'),title('Propeller revolutions n (rpm)'),grid

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)
