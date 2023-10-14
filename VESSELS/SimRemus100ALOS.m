% SIMremus100ALOS
% User editable script for 3-D path-following of the Remus 100 AUV 
% (remus100.m). The heading and pitch angles are controlled using two
% autopilots designed as PID controllers. The path is represented by 
% straigh lines going through waypoints. The heading and pitch commands
% are computed using the 3-D adaptive line-of-sight (ALOS) algorithm by
%
% T. I. Fossen and A. P. Aguiar (2023). A Uniform Semiglobal Exponential 
% Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path Following. 
% Automatica (submitted).
%
% Calls:      remus100.m     Remus 100 vehicle dynamics
%             ALOS3D.m       3-D adaptive LOS guidance law
%
% Author:     Thor I. Fossen
% Date:       2023-10-07
clearvars;
clear ALOS3D                 % clear the static variables used by ALOS3D

%% USER INPUTS
h  = 0.05;                  % sample time (s)
N  = 28000;                 % number of samples

% Ocean current speed and direction expressed in NED
Vc_0 = 0.5;                 % initial speed (m/s)
betaVc_0 = deg2rad(150);    % initial direction (rad)
w_c = 0.1;

% ALOS guidance law
Delta_h = 20;               % horizontal look-ahead distance (m)
Delta_v = 20;               % vertical look-ahead distance (m)
gamma_h = 0.002;            % adaptive gain, horizontal plane 
gamma_v = 0.002;            % adaptive gain, vertical plane 

M_theta = deg2rad(20);      % maximum value of estimates, alpha_c, beta_c
K_f = 0.3;                  % observer gain

% Waypoints
R_switch = 10;          % radius of sphere (m) used for waypoint switching 
wpt.pos.x = [0  -20  -100 0   200, 200  400];
wpt.pos.y = [0  200   600 950 1300 1800 2200];
wpt.pos.z = [0  10   100 100 50 50 50];

% Initial state vector: x = [ u v w p q r x y z eta eps1 eps2 eps3 ]'
phi = 0; theta = 0; psi = pi/2;     % Euler angles
U = 1.0;                            % surge velocity
quat = euler2q(phi,theta,psi);
x = [zeros(9,1); quat]; 
x(1) = U;

% Propeller initialization
n = 1000;               % initial propeller revolution (rpm)
n_d = 1300;             % desired propeller revolution, max 1525 rpm

% Autopilot integral states 
z_int = 0; theta_int = 0; psi_int = 0;

% Depth autopilot PID gains
Kp_theta = 2;            
Kd_theta = 3;
Ki_theta = 0.1;

% Heading autopilot PID gains
wn_b_psi = 0.9;             % bandwidth, pole placement algorithm 
wn_d_psi = wn_b_psi/5;      % desired natural frequency, reference model 
m66 = 7.5;                  % moment of inertia, yaw
Kp_psi = m66 * wn_b_psi^2;                 
Kd_psi = m66 * 2*wn_b_psi; 
Ki_psi = Kp_psi * (wn_b_psi/10);
   
%% MAIN LOOP
simdata = zeros(N+1,25); % allocate empty table for simulation data

for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % Measurements
   u = x(1);                % surge velocity (m/s)
   v = x(2);                % sway veloicty (m/s)
   q = x(5);                % pitch rate (rad/s)
   r = x(6);                % yaw rate (rad/s)
   xn = x(7);               % North position (m)
   yn = x(8);               % East position (m)   
   zn = x(9);               % Down position, depth (m)

   [phi,theta,psi] = q2euler(x(10:13));  % quaternion to Euler angles 

   % PID heading and depth controllers
   U_design = 2;        % ALOS design speed
   [psi_d, theta_d, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
       ALOS3D(xn,yn,zn,Delta_h,Delta_v,gamma_h,gamma_v,M_theta,h,...
       R_switch,wpt,U_design,K_f);       

   % Ocean current dynamics
   if t > 800 
       Vc_d = 0.65;
       w_V = 0.1;
       Vc = exp(-h*w_V) * Vc + (1 - exp(-h*w_V)) * Vc_d;
   else
       Vc = Vc_0;
   end

   if t > 500 
       betaVc_d = 160*pi/180;
       w_beta = 0.1;
       betaVc = exp(-h*w_beta) * betaVc + (1 - exp(-h*w_beta)) * betaVc_d;
   else
       betaVc = betaVc_0;
   end   
   
   betaVc = betaVc + (pi/180) * randn / 20;
   Vc = Vc + 0.002 * randn;
   
   % Autopilots for heading and pitch control
   delta_r = -Kp_psi * ssa( psi - psi_d ) - Kd_psi * r... 
       -Ki_psi * psi_int;                                          
   delta_s = -Kp_theta * ssa( theta - theta_d ) - Kd_theta * q ...
       -Ki_theta * theta_int;  

   % Propeller revolution (rpm)
   if (n < n_d)
       n = n + 1;
   end
   
   % Saturation, maximum controls
   max_ui = [30*pi/180 30*pi/180  1525]';   % rad, rad, rpm
   if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end  
   if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end  
   if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end
    
   ui = [delta_r delta_s n]';

   % Store simulation data in a table 
   simdata(i,:) = [t theta_d psi_d ui' x' y_e betaVc Vc beta_c_hat alpha_c_hat z_e];   

   % Propagate the vehicle dynamics (k+1)
   xdot = remus100(x,ui,Vc,betaVc,w_c);     % Remus 100 dyanmics
   x(1:9) = x(1:9) + h * xdot(1:9);         % Euler's integration method    
   quat = x(10:13);     
   quat = expm(Tquat(x(4:6)) * h) * quat;   % exact discretization
   x(10:13) = quat / norm(quat);            % normalization
   
   % Propagate integral states (k+1)
   theta_int = theta_int + h * ssa( theta - theta_d );
   psi_int = psi_int + h * ssa( psi - psi_d );   
   
end

%% PLOTS
t       = simdata(:,1);         
theta_d = simdata(:,2); 
psi_d   = simdata(:,3); 
u       = simdata(:,4:6); 
nu      = simdata(:,7:12);
y_e     = simdata(:,20);
betaVc  = simdata(:,21);
Vc      = simdata(:,22);
beta_c_hat = simdata(:,23);
alpha_c_hat = simdata(:,24);
z_e = simdata(:,25);

% Transform the unit quaternions to Euler angles
quaternion = simdata(:,16:19);
for i = 1:N+1
    [phi(i,1),theta(i,1),psi(i,1)] = q2euler(quaternion(i,:));   
end
eta = [simdata(:,13:15) phi theta psi];

alpha_c = atan( (nu(:,2).*sin(eta(:,4))+nu(:,3).*cos(eta(:,4))) ./ nu(:,1) );
Uv = nu(:,1) .* sqrt( 1 + tan(alpha_c).^2 );
beta_c = atan( ( nu(:,2).*cos(eta(:,4))-nu(:,3).*sin(eta(:,4)) ) ./ ... 
    ( Uv .* cos(eta(:,5)-alpha_c) ) );

uc = Vc .* cos(betaVc);
vc = Vc .* sin(betaVc);
wc = 0;
alpha = atan2( (nu(:,3)-wc), (nu(:,1)-uc) );
beta  = atan2( (nu(:,2)-vc), (nu(:,1)-uc) );
chi = eta(:,6) + beta(:,1);                     % course angle (rad)

%% Generalized velocity
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
set(findall(gcf,'type','legend'),'FontSize',12)

%% Control signals
figure(2); 
subplot(311),plot(t,(180/pi)*u(:,1))
xlabel('time (s)'),title('Rudder \delta_r (deg)'),grid
subplot(312),plot(t,(180/pi)*u(:,2))
xlabel('time (s)'),title('Stern planes \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3))
xlabel('time (s)'),title('Propeller revolutions n (rpm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)

%% Ocean current
figure(3); 
subplot(311),plot(t,nu(:,1),t,nu(:,2))
ylim([0,2.2])
xlabel('time (s)'),grid
legend('surge veloicty (m/s)','sway velocity (m/s)','Location','east')

subplot(312),plot(t,sqrt(nu(:,1).^2+nu(:,2).^2),t,Vc)
ylim([0,2.2])
xlabel('time (s)'),grid
legend('vehicle speed (m/s)','current speed (m/s)','Location','southeast')

subplot(313),plot(t,(180/pi)*betaVc)
xlabel('time (s)'),grid
legend('Current direction (deg)')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Autopilot performance
figure(4); 
subplot(211),plot(t,(180/pi)*eta(:,5),t,(180/pi)*theta_d)
xlabel('time (s)'),title('pitch angle (deg)'),grid
legend('true','desired')
subplot(212),plot(t,(180/pi)*eta(:,6),t,(180/pi)*psi_d)
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('true','desired')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)


%% Sideslip and angle of attack
figure(5); 

subplot(311)
plot(t,(180/pi)*alpha,'g',t,(180/pi)*alpha_c,'b',t,(180/pi)*alpha_c_hat,'r')
title('Angle of attack (deg)')
grid
legend('\alpha','\alpha_c','\alpha_c estimate')

subplot(312)
plot(t,(180/pi)*beta,'g',t,(180/pi)*beta_c,'b',t,(180/pi)*beta_c_hat,'r')
title('Sideslip angle (deg)')
xlabel('time (s)')
grid
legend('\beta','\beta_c','\beta_c estimate')

subplot(313)
plot(t,y_e,t,z_e)
title('Tracking errors (m)'),grid
xlabel('time (s)')
legend('cross-track error y_e^p','vertical-track error z_e^p')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Waypoints
figure(6); 
subplot(211)
plot(eta(:,2),eta(:,1))
hold on
plot(wpt.pos.y,wpt.pos.x,'rX','markersize',10)
hold off
xlabel('y^n (m)'),ylabel('x^n (m)')
xlim([0,2500])
axis('equal')
grid
legend('actual path','waypoints')

subplot(212)
plot(eta(:,2),eta(:,3))
hold on
plot(wpt.pos.y,wpt.pos.z,'rX','markersize',10)
hold off
xlim([0,2500])
ylim([0,150])
xlabel('y^n (m)'),ylabel('z^n (m)')
grid
legend('actual path','waypoints')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% 3-D position plot
figure(7); 
subplot(211)
plot(eta(:,2),eta(:,1))
title('North-East plot (m)')
xlabel('E'); ylabel('N'); grid
axis equal;
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

subplot(212)
plot3(eta(:,2),eta(:,1),eta(:,3))
hold on
plot3(wpt.pos.y,wpt.pos.x,wpt.pos.z,'ro','markersize',15)
hold off
title('North-East-Down plot (m)')
xlabel('E'); ylabel('N'); zlabel('D');
legend('actual path','waypoints'),grid
set(gca, 'ZDir', 'reverse');
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
view(-25, 30);  % view(AZ,EL) 

