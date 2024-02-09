% SIMremus100 User-editable script for simulating the Remus 100 AUV
% (remus100.m) under feedback control (simultaneously for depth and heading 
% control) when exposed to ocean currents. Several methods for heading and 
% depth control have been implemented, and switching between these methods 
% is done by specifying the control flags. Both the Euler angle and unit 
% quaternion representations of the Remus 100 model can be used. 
%
% Calls:      remus100.m (Remus 100 equations of motion)
%             integralSMCheading.m, refModel.m, q2euler.m, ssa.m
%
% Simulink:   demoAUVdepthHeadingControl.slx
%
% Author:     Thor I. Fossen
% Date:       2021-06-28
% Revisions:  2022-02-01 redesign of the autopilots
%             2022-05-06 retuning for new version of remus100.m
%             2022-05-08 added compability for unit quaternions
%             2023-02-09 added flags for multiple control laws

clearvars;
clear integralSMCheading    % reset integral state in integralSMCheading.m

%% USER INPUTS
h  = 0.02;                  % sample time (s)
N  = 25000;                 % number of samples

% Kinematic representation
kinematicsFlag = 1;         % 0: Euler angles
                            % 1: Unit quaternions
% Heading autopilot                   
headingControlFlag = 1;     % 0: PID pole placement algorithm
                            % 1: Intergral slidng mode control (SMC)

% Autopilot setpoints
n_d = 1525;                 % desired propeller revolution, max 1525 rpm
z_step = 30;                % step change in depth, max 100 m
psi_step = deg2rad(-60);    % step change in heading angle (rad)

mustBeInRange(n_d,0,1525);
mustBeInRange(z_step,0,100);

% Initial states
x = 0; y = 0; z = 10; phi = 0; theta = 0; psi = 0; 
n = 0;                      % initial propeller revolution (rpm)

% Autopilot integral states
z_int = 0;                  % depth
theta_int = 0;              % pitch angle
psi_int = 0;                % heading angle

% Intitial state vector
if (kinematicsFlag == 0)    % x = [ u v w p q r x y z phi theta psi ]'   
    x = [zeros(6,1); x; y; z; phi; theta; psi];
else                    % x = [ u v w p q r x y z eta eps1 eps2 eps3 ]'
    quat = euler2q(phi,theta,psi);
    x = [zeros(6,1); x; y; z; quat];
end

% Ocean current speed and direction expressed in NED
Vc = 0.5;                   % speed (m/s)
betaVc = deg2rad(170);      % direction (rad)

% Depth controller (suceessive-loop closure)
z_d = z;                    % initial depth (m), reference model
wn_d_z = 0.02;              % desired natural frequency, reference model
Kp_z = 0.1;                 % proportional gain (heave)
T_z = 100;                  % integral time constant (heave)
Kp_theta = 5.0;             % proportional gain (pitch)
Kd_theta = 2.0;             % derivative gain (pitch)
Ki_theta = 0.3;             % integral gain (pitch
K_w = 5.0;                  % optional heave velocity feedback gain

% Nomoto gain parameters
K_yaw = 5/20;               % K_yaw = r_max / delta_max
T_yaw = 1;                  % Time constant in yaw

% Heading autopilot reference model 
psi_d = psi;                % initial heading angle (rad), reference model
r_d = 0;                    % initial yaw rate (rad/s), reference model
a_d = 0;                    % initial yaw acc. (rad/s^2), reference model
zeta_d_psi = 1.0;           % desired relative damping factor, refence model
wn_d_psi = 0.1;             % desired natural frequency, reference model 
r_max = deg2rad(5.0);       % maximum turning rate (rad/s)

% Heading autopilot (Equation 16.479 in Fossen 2021)
% sigma = r-r_d + 2*lambda*ssa(psi-psi_d) + lambda^2 * integral(ssa(psi-psi_d))
% delta = (T_yaw*r_r_dot + r_r - K_d*sigma - K_sigma*(sigma/phi_b)) / K_yaw
lambda = 0.1;
phi_b = 0.1;                % boundary layer thickness

if headingControlFlag == 0  % PID controller
    K_d = 0.5; 
    K_sigma = 0; 
else                        % SMC controller 
    K_d = 0;
    K_sigma = 0.05;             
end
   
%% Display
disp('-------------------------------------------------------------');
disp('MSS toolbox: Remus 100 AUV (Length = 1.6 m, Diameter = 19 cm)')  
if (kinematicsFlag == 0)
    disp('Euler angle representation (12 states)')
else
    disp('Unit quaternion representation (13 states)') 
end
if (headingControlFlag == 0)
    disp('Heading autopilot: PID poleplacement control')
else
    disp('Heading autopilot: Intergral slidng mode control (SMC) ')   
end
disp('Depth autopilot:   Succesive-loop closure')
disp('-------------------------------------------------------------');

%% MAIN LOOP
simdata = zeros(N+1,length(x)+8); % allocate empty table for simulation data

for i = 1:N+1
    
   t = (i-1)*h;                 % time
   
   % Measurements
   w = x(3);                    % heave velocity
   q = x(5);                    % pitch rate
   r = x(6);                     % yaw rate
   z = x(9);                    % z-position (depth)
   Uv = sqrt(x(1)^2+x(3)^2);    % vertical speed
   
   if (kinematicsFlag==0)
         phi = x(10); theta = x(11); psi = x(12);   % Euler angles
   else
         [phi,theta,psi] = q2euler(x(10:13));  % quaternion to Euler angles
   end
   
   % Depth command, z_ref  
   if t > 200
       z_ref = z_step;
   else
       z_ref = 10;
   end

   % LP filtering the depth command
   if Uv < 1.0           % reduce bandwidth at low speed
       wnz = wn_d_z / 2; 
   else
       wnz = wn_d_z;
   end  
   z_d = exp(-h*wnz) * z_d + (1 - exp(-h*wnz)) * z_ref; 

   % Depth autopilot using the stern planes (succesive-loop closure)
   theta_d = Kp_z * ( (z - z_d) + (1/T_z) * z_int );            % PI 
   delta_s = -Kp_theta * ssa( theta - theta_d )...              % PID
        - Kd_theta * q - Ki_theta * theta_int - K_w * w;                                        

   % PID heading angle command, psi_ref
   if t > 200
      psi_ref = psi_step;
   else
      psi_ref = deg2rad(0);       
   end   
   
   % Third-order reference model for the heading angle   
   [psi_d,r_d,a_d] = refModel(psi_d,r_d,a_d,psi_ref,r_max,...
        zeta_d_psi,wn_d_psi,h,1);
  
   % Heading autopilot using the tail rudder
    delta_r = integralSMCheading(psi,r,psi_d,r_d,a_d,K_d,K_sigma,1,...
          phi_b,K_yaw,T_yaw,h);
   
   % Propeller revolution (rpm)
   if (n < n_d)
       n = n + 1;
   end

    % Amplitude saturation of the control signals
    n_max = 1525;                                % maximum propeller rpm
    max_ui = [deg2rad(15) deg2rad(15) n_max]';   % deg, deg, rpm

    if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
    if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
    if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end

    ui = [delta_r delta_s n]';                % Commanded control inputs 
   
   % Store simulation data in a table 
   simdata(i,:) = [t z_d theta_d psi_d r_d ui' x'];   
   
   % Propagate the vehicle dynamics (k+1)
   xdot = remus100(x,ui,Vc,betaVc);
   
   if (kinematicsFlag == 0)     
       x = x + h * xdot;                       % Euler's integration method
   else
       x(1:9) = x(1:9) + h * xdot(1:9);        % Euler's integration method       
       quat = x(10:13);                        % unit quaternion
       quat = expm(Tquat(x(4:6)) * h) * quat;  % exact discretization
       x(10:13) = quat / norm(quat);           % normalization
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
r_d     = simdata(:,5); 
u       = simdata(:,6:8); 
nu      = simdata(:,9:14);

if (kinematicsFlag==0)         % Euler angle representation
    eta = simdata(:,15:20);
else                       % Transform the unit quaternions to Euler angles
    quaternion = simdata(:,18:21);
    for i = 1:N+1
        [phi(i,1),theta(i,1),psi(i,1)] = q2euler(quaternion(i,:));
    end
    eta = [simdata(:,15:17) phi theta psi];
    
    figure(4)
    plot(t,quaternion,'linewidth',2);
    xlabel('time (s)'),title('Unit quaternion'),grid
    legend('\eta','\epsilon_1','\epsilon_2','\epsilon_3')
end

%% Generalized velocity
figure(1); 
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
legend('true','desired')
xlabel('time (s)'),title('Yaw rate (deg/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

%% Speed, heave position and Euler angles
figure(2); 
subplot(511),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2+nu(:,3).^2));
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(512),plot(t,eta(:,3),t,z_d)
xlabel('time (s)'),title('heave position (m)'),grid
legend('true','desired')
subplot(513),plot(t,rad2deg(eta(:,4)))
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(514),plot(t,rad2deg(eta(:,5)),t,rad2deg(theta_d))
xlabel('time (s)'),title('pitch angle (deg)'),grid
legend('true','desired')
subplot(515),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d))
xlabel('time (s)'),title('yaw angle (deg)'),grid
legend('true','desired')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

%% Control signals
figure(3); 
subplot(311),plot(t,rad2deg(u(:,1)))
xlabel('time (s)'),title('Rudder command \delta_r (deg)'),grid
subplot(312),plot(t,rad2deg(u(:,2)))
xlabel('time (s)'),title('Stern-plane command \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3))
xlabel('time (s)'),title('Propeller speed command n (rpm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

%% 3-D position plot
figure(4); 
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
title('North-East-Down plot (m)')
xlabel('E'); ylabel('N'); zlabel('D'); grid
set(gca, 'ZDir', 'reverse');
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
view(-25, 30);  % view(AZ,EL) 


