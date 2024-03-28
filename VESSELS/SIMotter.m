% SIMotter User editable script for simulation of the Otter USV ('otter.m')
% under feedback control when exposed to ocean currents. Methods for
% heading and path-following control have been implemented, and switching 
% between these methods is done by specifying the control flag. 
%
% Calls:      otter.m
%             refModel.m
%             ALOS.m
%
% Author:    Thor I. Fossen
% Date:      2021-04-25
% Revisions: 2023-10-14 Added a heading autopilot and reference model
%            2024-03-28 Added ALOS path-following control

clear ALOSpsi
clearvars;

%% USER INPUTS
h  = 0.02;        % sampling time [s]
N  = 6000;		  % number of samples

% Control system                  
ControlFlag = 1;                % 0: PID heading autopilot
                                % 1: ALOS path-following control   

% Load condition
mp = 25;                        % payload mass (kg), max value 45 kg
rp = [0.05 0 -0.35]';           % location of payload (m)

% Ocean current
V_c = 0.3;                      % current speed (m/s)
beta_c = deg2rad(30);           % current direction (rad)

% Waypoint path-following parameters
U = 1.0;                        % cruise speed
Delta = 10;                     % look-ahead distance
gamma = 0.001;                  % adaptive gain
R_switch = 3;                   % radius of switching circle
K_f = 0.5;                      % desired yaw rate observer gain
wpt.pos.x = [0,50,80,80]';      % waypoints
wpt.pos.y = [0,50,100,200]';

% PID heading autopilot (Nomoto gains)
T = 1;
m = 41.4;          % m = T/K
K = T / m;

wn = 2.5;          % PID pole-placement parameters
zeta = 1.0;

Kp = m * wn^2;
Kd = m * (2 * zeta * wn - 1/T);
Td = Kd / Kp; 
Ti = 10 / wn;

% Reference model parameters
wn_d = 1.0;                     % natural frequency (rad/s)
zeta_d = 1.0;                   % relative damping factor (-)
r_max = deg2rad(5.0);           % maximum turning rate (rad/s)

% Input matrix
y_prop = 0.395;                 % distance from centerline to propeller (m)
k_pos = 0.0111;                 % positive Bollard, one propeller 
B = k_pos * [ 1 1               % input matrix
             y_prop  -y_prop  ];
Binv = inv(B);

% Propeller dynamics
T_n = 0.1;                      % propeller time constant (s)      
n = [0 0]';                     % n = [ n_left n_right ]'

% Initial states
eta = [0 0 0 0 0 0]';           % eta = [x y z phi theta psi]' 
nu  = [0 0 0 0 0 0]';           % nu  = [u v w p q r]'	 
z_psi = 0;                      % integral state
psi_d = eta(6);                 % reference model states
r_d = 0;
a_d = 0;

%% Display
disp('-------------------------------------------------------------');
disp('MSS toolbox: Otter USV (Length = 1.6 m, Beam = 1.08 m)')  
if (ControlFlag == 0)
    disp('PID heading autopilot with reference feeforward')
else
    disp('ALOS guidance law for path-following control')   
end
disp('-------------------------------------------------------------');

%% MAIN LOOP
simdata = zeros(N+1,15);                   % table for simulation data

for i=1:N+1

   t = (i-1) * h;                          % time (s)    
   
   % Measurements
   x = eta(1) + 0.001 * randn;
   y = eta(2) + 0.001 * randn;
   r = nu(6) + 0.001 * randn;
   psi = eta(6) + 0.001 * randn;
 
   % Heading commands
   if ControlFlag == 0 % Heading autopilot
       if ( t < 20 )
           psi_ref = deg2rad(20);
       else
           psi_ref = deg2rad(0);
       end
   else % path following
       [psi_d,r_d] = ALOSpsi(x,y,Delta,gamma,h,R_switch,wpt,U,K_f);
       a_d = 0;
   end

   % Heading autopilot
   tau_X = 100;
   tau_N = (T/K) * a_d + (1/K) * r_d -...
       Kp * (ssa( psi-psi_d) +...
       Td * (r - r_d) + (1/Ti) * z_psi );
   u = Binv * [tau_X tau_N]';
   n_c = sign(u) .* sqrt( abs(u) );

   % Store simulation data in a table   
   simdata(i,:) = [t eta' nu' r_d psi_d];    
   
   % USV dynamics
   xdot = otter([nu; eta],n,mp,rp,V_c,beta_c);

   % Euler's integration method (k+1)
   nu = nu + h * xdot(1:6);               
   eta = eta + h * xdot(7:12); 
   n = n + h/T_n * (n_c - n);              
   z_psi = z_psi + h * ssa( psi-psi_d );  
   
   % Propagation of reference model(k+1)
   if ControlFlag == 0
       [psi_d,r_d,a_d] = refModel(psi_d,r_d,a_d,psi_ref,r_max,zeta_d,wn_d,h,1);                  
   end

end

%% PLOTS
t = simdata(:,1); 
eta = simdata(:,2:7); 
nu  = simdata(:,8:13); 
r_d = simdata(:,14); 
psi_d = simdata(:,15); 

clf

figure(1); figure(gcf)
plot(eta(:,2),eta(:,1)); 
if ControlFlag == 1
    hold on; plot(wpt.pos.y,wpt.pos.x, 'rx'); hold off
end
xlabel('East (m)', 'FontSize', 14);
ylabel('North (m)', 'FontSize', 14);
title('North-East positions', 'FontSize', 14);
axis equal; grid
set(findall(gcf,'type','line'),'linewidth',2)

figure(2); figure(gcf)
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
legend('r','r_d','Location','best')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(3); figure(gcf)
subplot(511),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2));
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(512),plot(t,eta(:,3),'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
subplot(513),plot(t,rad2deg(eta(:,4)))
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(514),plot(t,rad2deg(eta(:,5)))
xlabel('time (s)'),title('pitch angle (deg)'),grid
subplot(515),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d))
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('\psi','\psi_d','Location','best')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
