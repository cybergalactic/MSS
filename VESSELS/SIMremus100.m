% SIMremus100 
% User-editable script for simulating the Remus 100 AUV ('remus100.m') under
% depth and heading control when exposed to ocean currents. The autpilots 
% can also be combined with the 3-D Adaptive Line-Of-Sight (ALOS) guidance 
% law by Fossen and Aguiar (2024) for path following. The ALOS algorithm
% computes the desired heading and pitch angles when the path is a straight
% line segment going through the waypoints (wpt.pos.x,wpt.pos.y,wpt.pos.z).
% Both the Euler angle and unit quaternion representations of the Remus 100 
% model can be used. The following methods are available for selection:
%
% 1: Heading control: PID pole-placement algorithm
% 2. Heading control: Integral slidng mode control (SMC)
% 3. Path-following:  ALOS guidance law for 3-D path following
%
% Calls:    remus100.m (Remus 100 equations of motion)
%           integralSMCheading.m
%           refModel.m
%           ALOS.m
%           LOSobserver.m
%           q2euler.m
%           ssa.m
%
% Simulink: demoAUVdepthHeadingControl.slx
%
% Reference:
%   T. I. Fossen and P. Aguiar (2024). A Uniform Semiglobal Exponential 
%   Stable Adaptive Line-of-Sight (ALOS) Guidance Law for 3-D Path Following. 
%   Automatica 163, 111556. doi.org/10.1016/j.automatica.2024.111556
%
% Author:     Thor I. Fossen
% Date:       2021-06-28
% Revisions:  2022-02-01 Redesign of the autopilots
%             2022-05-06 Retuning for new version of remus100.m
%             2022-05-08 Added compability for unit quaternions
%             2024-03-27 Using forward and backward Euler to integrate xdot
%             2024-04-02 Added simulation menu and ALOS path following

clearvars;                  % clear variables from memory
clear integralSMCheading    % clear the integral state in integralSMCheading
clear ALOS3D                % clear the static variables used by ALOS3D
close all;                  % closes all figure windows 

%% USER INPUTS
h  = 0.05;                  % sample time (s)
N  = 28000;                 % number of samples

[ControlFlag, KinematicsFlag] = controlMethod(); % choose control method

% Waypoints
wpt.pos.x = [0  -20 -100   0  200, 200  400];
wpt.pos.y = [0  200  600 950 1300 1800 2200];
wpt.pos.z = [0   10  100 100   50   50   50];

% Initial states
xn = 0; yn = 0; zn = 0;         % initial NED positions
phi = 0;                        % intial Euler angles
theta = 0; 
psi = atan2(wpt.pos.y(2)-wpt.pos.y(1),wpt.pos.x(2)-wpt.pos.x(1)); 
U = 1;                          % initial speed 
theta_d = 0; q_d = 0;           % initial pitch reference signals
psi_d = psi; r_d = 0; a_d = 0;  % initial yaw reference signals


% Intitial state vector
if (KinematicsFlag == 1)  % x = [ u v w p q r x y z phi theta psi ]'   
    x = [U; zeros(5,1); xn; yn; zn; phi; theta; psi];
else                      % x = [ u v w p q r x y z eta eps1 eps2 eps3 ]'
    quat = euler2q(phi,theta,psi);
    x = [U; zeros(5,1); xn; yn; zn; quat];
end

% Ocean current speed and direction expressed in NED
Vc = 0.5;                   % horisontal speed (m/s)
betaVc = deg2rad(30);       % horizontal direction (rad)
wc = 0.1;                   % vertical speed (m/s)

% Propeller initialization
n = 1000;                   % initial propeller revolution (rpm)
n_d = 1300;                 % desired propeller revolution, max 1525 rpm
mustBeInRange(n_d,0,1525);

%% HEADING AND DEPTH AUTOPILOTS PARAMETERS
psi_step = deg2rad(-60);    % step change in heading angle (rad)
z_step = 30;                % step change in depth, max 100 m
mustBeInRange(z_step,0,100);

% Autopilot integral states
z_int = 0;                  % depth
theta_int = 0;              % pitch angle
psi_int = 0;                % heading angle

% Depth controller (suceessive-loop closure)
z_d = zn;                   % initial depth (m), reference model
wn_d_z = 0.02;              % desired natural frequency, reference model
Kp_z = 0.1;                 % proportional gain (heave)
T_z = 100;                  % integral time constant (heave)
Kp_theta = 5.0;             % proportional gain (pitch)
Kd_theta = 2.0;             % derivative gain (pitch)
Ki_theta = 0.3;             % integral gain (pitch
K_w = 5.0;                  % optional heave velocity feedback gain

% Feedforward gains (Nomoto gain parameters)
K_yaw = 5/20;               % K_yaw = r_max / delta_max
T_yaw = 1;                  % Time constant in yaw

% Heading autopilot reference model 
zeta_d_psi = 1.0;           % desired relative damping factor, refence model
wn_d_psi = 0.1;             % desired natural frequency, reference model 
r_max = deg2rad(5.0);       % maximum turning rate (rad/s)

% Heading autopilot (Equation 16.479 in Fossen 2021)
% sigma = r-r_d + 2*lambda*ssa(psi-psi_d) + lambda^2 * integral(ssa(psi-psi_d))
% delta = (T_yaw*r_r_dot + r_r - K_d*sigma - K_sigma*(sigma/phi_b)) / K_yaw
lambda = 0.1;
phi_b = 0.1;                % boundary layer thickness

if ControlFlag == 1         % PID controller
    K_d = 0.5; 
    K_sigma = 0; 
else                        % SMC controller 
    K_d = 0;
    K_sigma = 0.05;             
end

%% ALOS PATH-FOLLOWING PARAMETERS
Delta_h = 20;               % horizontal look-ahead distance (m)
Delta_v = 20;               % vertical look-ahead distance (m)
gamma_h = 0.001;            % adaptive gain, horizontal plane 
gamma_v = 0.001;            % adaptive gain, vertical plane 
M_theta = deg2rad(20);      % maximum value of estimates, alpha_c, beta_c

% Additional parameter for straigh-line path following
R_switch = 5;               % radius of switching circle
K_f = 0.5;                  % pitch and yaw rate observer gain

%% MAIN LOOP
simdata = zeros(N+1,length(x)+11); % allocate empty table for simulation data
ALOSdata = zeros(N+1,4);           % allocate empty table for ALOS data

for i = 1:N+1
    
   t = (i-1)*h;             % time
   
   % Measurements
   u = x(1);                % surge velocity (m/s)
   v = x(2);                % sway veloicty (m/s)
   w = x(3);                % heave veloicty (m/s)   
   q = x(5);                % pitch rate (rad/s)
   r = x(6);                % yaw rate (rad/s)
   xn = x(7);               % North position (m)
   yn = x(8);               % East position (m)   
   zn = x(9);               % Down position, depth (m)
   
   if (KinematicsFlag == 1)
         phi = x(10); theta = x(11); psi = x(12); % Euler angles
   else
         [phi,theta,psi] = q2euler(x(10:13)); % quaternion to Euler angles
   end
   
   % Control systems
   if ControlFlag == 1 || ControlFlag == 2

       % Depth command, z_ref
       if t > 200
           z_ref = z_step;
       else
           z_ref = 10;
       end

       % LP filtering the depth command
       Uv = sqrt(x(1)^2+x(3)^2);    % vertical speed
       if Uv < 1.0                  % reduce bandwidth at low speed
           wnz = wn_d_z / 2;
       else
           wnz = wn_d_z;
       end
       z_d = exp(-h*wnz) * z_d + (1 - exp(-h*wnz)) * z_ref;

       % Depth autopilot using the stern planes (succesive-loop closure)
       theta_d = Kp_z * ( (zn - z_d) + (1/T_z) * z_int );     % PI
       delta_s = -Kp_theta * ssa( theta - theta_d )...        % PID
           - Kd_theta * q - Ki_theta * theta_int - K_w * w;

       % PID heading angle command, psi_ref
       if t > 200
           psi_ref = psi_step;
       else
           psi_ref = deg2rad(0);
       end

       % Heading autopilot using the tail rudder
       delta_r = integralSMCheading(psi, r, psi_d, r_d, a_d, ...
           K_d, K_sigma, 1, phi_b, K_yaw, T_yaw, h);

       % Third-order reference model for the heading angle
       [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
           zeta_d_psi, wn_d_psi, h, 1);

   else % ControlFlag = 3

       % Heading autopilot using the tail rudder (integral SMC)
       delta_r = integralSMCheading(psi, r, psi_d, r_d, a_d, ...
           K_d, K_sigma, 1, phi_b, K_yaw, T_yaw, h);
       
       % Depth autopilot using the stern planes (PID)
       delta_s = -Kp_theta * ssa( theta - theta_d )...             
           - Kd_theta * q - Ki_theta * theta_int - K_w * w;

       % ALOS guidance law
       [psi_ref, theta_ref, y_e, z_e, alpha_c_hat, beta_c_hat] = ...
           ALOS3D(xn, yn, zn, Delta_h, Delta_v, gamma_h, gamma_v,...
           M_theta, h, R_switch, wpt);

       % ALOS observer
       [theta_d, q_d] = LOSobserver(theta_d, q_d, theta_ref, h, K_f);
       [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);
       if abs(r_d) > r_max, r_d = sign(r_d) * r_max; end

       % Ocean current dynamics
       if t > 800
           Vc_d = 0.65;
           w_V = 0.1;
           Vc = exp(-h*w_V) * Vc + (1 - exp(-h*w_V)) * Vc_d;
       else
           Vc = 0.5;
       end

       if t > 500
           betaVc_d = deg2rad(160);
           w_beta = 0.1;
           betaVc = exp(-h*w_beta) * betaVc + (1 - exp(-h*w_beta)) * betaVc_d;
       else
           betaVc = deg2rad(150);
       end

       betaVc = betaVc + (pi/180) * randn / 20;
       Vc = Vc + 0.002 * randn;

       ALOSdata(i,:) = [y_e z_e alpha_c_hat beta_c_hat];

   end

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
   simdata(i,:) = [t z_d theta_d psi_d r_d Vc betaVc wc ui' x'];   
   
   % Propagate the vehicle dynamics (k+1), (Fossen 2021, Eq. B27-B28)
   % x = x + h * xdot is replaced by forward and backward Euler integration
   xdot = remus100(x, ui, Vc, betaVc, wc);

   if (KinematicsFlag == 1)     
       Jmtrx = eulerang(x(10),x(11),x(12));
       x(1:6) = x(1:6) + h * xdot(1:6);        % forward Euler 
       x(7:12) = x(7:12) + h * Jmtrx * x(1:6); % backward Euler 
   else
       quat = x(10:13);                        % unit quaternion
       Rq = Rquat(quat);                       % rotation matrix
       x(1:6) = x(1:6) + h * xdot(1:6);        % forward Euler       
       x(7:9) = x(7:9) + h * Rq * x(1:3);      % backward Euler
       quat = expm(Tquat(x(4:6)) * h) * quat;  % exact quat. discretization
       x(10:13) = quat / norm(quat);           % normalization
   end
   
   % Euler's integration method (k+1)
   z_int = z_int + h * ( zn - z_d );
   theta_int = theta_int + h * ssa( theta - theta_d );
   psi_int = psi_int + h * ssa( psi - psi_d );   
   
end

%% PLOTS
screenSize = get(0, 'ScreenSize'); % Returns [left bottom width height]
screenW = screenSize(3); 
screenH = screenSize(4);

% simdata = [t z_d theta_d psi_d r_d Vc betaVc wc ui' x']
t       = simdata(:,1);  
z_d     = simdata(:,2);
theta_d = simdata(:,3); 
psi_d   = simdata(:,4); 
r_d     = simdata(:,5); 
Vc      = simdata(:,6);
betaVc  = simdata(:,7);
wc      = simdata(:,8);
u       = simdata(:,9:11); 
nu      = simdata(:,12:17);

if (KinematicsFlag == 1) % Euler angle representation
    eta = simdata(:,18:23);
else % Transform the unit quaternions to Euler angles
    quaternion = simdata(:,21:24);
    for i = 1:N+1
        [phi(i,1),theta(i,1),psi(i,1)] = q2euler(quaternion(i,:));
    end
    eta = [simdata(:,18:20) phi theta psi];
end

% ALOSdata = [y_e z_e alpha_c_hat beta_c_hat]
y_e = ALOSdata(:,1);
z_e = ALOSdata(:,2);
alpha_c_hat = ALOSdata(:,3);
beta_c_hat = ALOSdata(:,4);

uc = Vc .* cos(betaVc);
vc = Vc .* sin(betaVc);
alpha_c = atan( (nu(:,2).*sin(eta(:,4))+nu(:,3).*cos(eta(:,4))) ./ nu(:,1) );
Uv = nu(:,1) .* sqrt( 1 + tan(alpha_c).^2 );
beta_c = atan( ( nu(:,2).*cos(eta(:,4))-nu(:,3).*sin(eta(:,4)) ) ./ ... 
    ( Uv .* cos(eta(:,5)-alpha_c) ) );
alpha = atan2( (nu(:,3)-wc), (nu(:,1)-uc) );
beta  = atan2( (nu(:,2)-vc), (nu(:,1)-uc) );

%% Generalized velocity
figure(1); set(gcf, 'Position', [1, 1, screenW/3, screenH]); 
subplot(611),plot(t,nu(:,1))
xlabel('Time (s)'),title('Surge velocity (m/s)'),grid
subplot(612),plot(t,nu(:,2))
xlabel('Time (s)'),title('Sway velocity (m/s)'),grid
subplot(613),plot(t,nu(:,3))
xlabel('Time (s)'),title('Heave velocity (m/s)'),grid
subplot(614),plot(t,(180/pi)*nu(:,4))
xlabel('Time (s)'),title('Roll rate (deg/s)'),grid
subplot(615),plot(t,(180/pi)*nu(:,5))
xlabel('Time (s)'),title('Pitch rate (deg/s)'),grid
subplot(616),plot(t,(180/pi)*nu(:,6))
xlabel('Time (s)'),title('Yaw rate (deg/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)

%% Heave position and Euler angles
figure(2); set(gcf, 'Position', [screenW/3, 1, screenW/3, screenH]);
if ControlFlag == 3; z_d = eta(:,3); end
subplot(411),plot(t,eta(:,3),t,z_d)
xlabel('Time (s)'),title('Heave position (m)'),grid
legend('True','Desired')
subplot(412),plot(t,rad2deg(eta(:,4)))
xlabel('Time (s)'),title('Roll angle (deg)'),grid
subplot(413),plot(t,rad2deg(eta(:,5)),t,rad2deg(theta_d))
xlabel('Time (s)'),title('Pitch angle (deg)'),grid
legend('True','Desired')
subplot(414),plot(t,rad2deg(eta(:,6)),t,rad2deg(psi_d))
xlabel('Time (s)'),title('Yaw angle (deg)'),grid
legend('True','Desired')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

%% Control signals
figure(3); set(gcf,'Position',[2*screenW/3,screenH/2,screenW/3,screenH/2]);
subplot(311),plot(t,rad2deg(u(:,1)))
xlabel('Time (s)'),title('Rudder command \delta_r (deg)'),grid
subplot(312),plot(t,rad2deg(u(:,2)))
xlabel('Time (s)'),title('Stern-plane command \delta_s (deg)'),grid
subplot(313),plot(t,u(:,3))
xlabel('Time (s)'),title('Propeller speed command n (rpm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',16)

%% Ocean currents and speed
figure(4); set(gcf,'Position',[2*screenW/3,1,screenW/3,screenH/2]);
subplot(311),plot(t,sqrt(nu(:,1).^2+nu(:,2).^2),t,Vc)
xlabel('Time (s)'),grid
legend('Vehicle horizontal speed (m/s)','Ocean current horizontal speed (m/s)',...
    'Location','best')
subplot(312),plot(t,nu(:,3),t,wc)
xlabel('Time (s)'),grid
legend('Vehicle heave velocity (m/s)','Ocean current heave velcoity (m/s)',...
    'Location','best')
subplot(313),plot(t,rad2deg(betaVc),'r')
xlabel('Time (s)'),grid
legend('Ocean current direction (deg)','Location','best')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Sideslip and angle of attack
if ControlFlag == 3
    figure(5); set(gcf,'Position',[100,100,screenW/3,screenH]);
    subplot(311)
    plot(t,rad2deg(alpha),'g',t,rad2deg(alpha_c),'b',...
        t,rad2deg(alpha_c_hat),'r')
    title('Angle of attack (deg)')
    xlabel('Time (s)')
    grid
    legend('\alpha','\alpha_c','\alpha_c estimate','Location','best')
    subplot(312)
    plot(t,rad2deg(beta),'g',t,rad2deg(beta_c),'b',...
        t,rad2deg(beta_c_hat),'r')
    title('Sideslip angle (deg)')
    xlabel('Time (s)')
    grid
    legend('\beta','\beta_c','\beta_c estimate','Location','best')
    subplot(313)
    plot(t,y_e,t,z_e)
    title('Tracking errors (m)'),grid
    xlabel('Time (s)')
    legend('Cross-track error y_e^p','Vertical-track error z_e^p')
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
end

%% 2-D position plots with waypoints
if ControlFlag == 3
    figure(6); set(gcf,'Position',[300,200,screenW/3,screenH/2]);
    subplot(211)
    plot(eta(:,2),eta(:,1))
    hold on; 
    plot(wpt.pos.y,wpt.pos.x,'rX','markersize',10); 
    hold off
    xlabel('East'), ylabel('North')
    title('North-East plot (m)')
    xlim([0,2500])
    axis('equal')
    grid
    legend('Actual path','Waypoints','Location','best')
    subplot(212)
    plot(eta(:,2),eta(:,3))
    hold on; 
    plot(wpt.pos.y,wpt.pos.z,'rX','markersize',10); 
    hold off
    xlim([0,2500])
    ylim([0,150])
    xlabel('East'),ylabel('Down')
    title('Down-East plot (m)')
    grid
    legend('Actual path','Waypoints','Location','best')
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
end

%% 3-D position plot with waypoints
if ControlFlag == 3
    figure(7); 
    plot3(eta(:,2),eta(:,1),eta(:,3))
    hold on; 
    plot3(wpt.pos.y,wpt.pos.x,wpt.pos.z,'ro','markersize',15); 
    hold off
    title('North-East-Down plot (m)')
    xlabel('East'); ylabel('North'); zlabel('Down');
    legend('Actual path','Waypoints','Location','best'),grid
    set(gca, 'ZDir', 'reverse');
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    view(-25, 30);  % view(AZ,EL)
end


%% DISPLAY AND CHOOSE CONTROL METHOD
function [ControlFlag, KinematicsFlag] = controlMethod()

    ControlFlag = 1;        % Default to 1 for "PID pole-placement control"
    KinematicsFlag = 1;     % Default to 1 for "Euler angles"

    fontSize = 13;

    f = figure('Position', [400, 400, 570, 300], 'Name', ...
        'Select Control Method and Kinematic Representation', ...
        'MenuBar','none','NumberTitle','off','WindowStyle','modal');

    uicontrol('Style', 'text', 'String','Control Method:', ...
        'Position', [10 280 180 15],'HorizontalAlignment','left',...
        'FontSize',fontSize,'FontWeight','bold');

    % Add button group for control methods
    bg = uibuttongroup('Parent', f, 'Position', [0.02 0.6 0.96 0.3], ...
        'SelectionChangedFcn', @controlSelection);
    
    % Add radio buttons within the button group
    uicontrol(bg, 'Style', 'radiobutton', 'FontSize', fontSize, 'String', ...
        'Heading autopilot: PID pole-placement control', ...
        'Position', [10 55 480 30], 'HandleVisibility', 'on', 'Tag', '1');
    uicontrol(bg, 'Style', 'radiobutton', 'FontSize', fontSize, 'String', ...
        'Heading autopilot: Integral sliding mode control (SMC)', ...
        'Position', [10 30 480 30], 'HandleVisibility', 'on', 'Tag', '2');
    uicontrol(bg, 'Style', 'radiobutton','FontSize', fontSize, 'String', ...
        ['3-D path-following: Adaptive Line-Of-Sight (ALOS) guidance ' ...
        'law for heading control'], 'Position', [10 5 560 30], ...
        'HandleVisibility', 'on', 'Tag', '3');

    uicontrol('Style','text','FontSize',fontSize,'FontWeight','bold','String',...
        'Choose Euler Angles or Unit Quaternions to Run the Simulation:', ...
        'Position',[10 140 480 20],'HorizontalAlignment','left');

    % Add toggle buttons for Kinematics options
    b1 = uicontrol('Style', 'togglebutton','FontSize', fontSize,...
        'String', 'Euler angles','Position', [50 100 150 30], ...
        'Callback', @toggleButtonCallback,'UserData',1,'Tag','1');
    b2 = uicontrol('Style', 'togglebutton','FontSize', fontSize, ...
        'String', 'Unit quaternions','Position', [50 65 150 30], ...
        'Callback', @toggleButtonCallback,'UserData',2,'Tag','2');

    uiwait(f); % Wait for the user to make a selection and close the figure

    % Callback function for toggle buttons
    function toggleButtonCallback(src, ~)
        % Ensure only one toggle button can be selected at a time
        set([b1, b2], 'Value', 0); % Reset both buttons
        src.Value = 1; % Set the clicked one to selected
        KinematicsFlag = src.UserData; % Update KinematicsFlag based on selection
        runSimulation(); % Run the simulation immediately after selection
    end

    % Selection change function for radio buttons
    function controlSelection(~, event)
        ControlFlag = str2double(event.NewValue.Tag);
    end

    % Function to run simulation and display choice
    function runSimulation()
        uiresume(f); % Resume execution (to return values) before closing
        close(f); 
        
        % Display the selected configuration
        disp('-------------------------------------------------------------');
        disp('MSS toolbox: Remus 100 AUV (Length = 1.6 m, Diameter = 19 cm)');
        if (KinematicsFlag == 1)
            disp('Euler angle representation (12 states)');
        else
            disp('Unit quaternion representation (13 states)');
        end
        disp(' ')
        if (ControlFlag == 1)
            disp('Heading autopilot: PID poleplacement control');
            disp('Depth autopilot:   Successive-loop closure');
        elseif (ControlFlag == 2)
            disp('Heading autopilot: Integral sliding mode control (SMC)');
            disp('Depth autopilot:   Successive-loop closure');
        else
            disp('Heading autopilot: Integral sliding mode control (SMC)');
            disp('Depth autopilot:   Successive-loop closure');
            disp('Path-following:    ALOS guidance law for 3-D path following')
        end
        disp('-------------------------------------------------------------');
    end
end