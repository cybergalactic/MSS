% SIMotter: User editable script for simulation of the Otter USV('otter.m'),
% length 2.0 m, beam 1.08 m, and mass 55 kg under feedback control when 
% exposed to ocean currents; see 'otter.m'. Methods for heading and 2-D
% adaptive line-of-sight (ALOS) and integral line-of-sight (ILOS) 
% path-following control have been implemented, and switching between 
% these methods is done by specifying the flag:
%
% ControlFlag: 0 - PID heading autopilot, no path following
%              1 - ALOS path-following control using straight lines 
%                  and waypoint switching
%              2 - ILOS path-following control using straight lines 
%                  and waypoint switching
%              3 - ALOS path-following control using Hermite splines
%
% Calls:      otter.m
%             refModel.m
%             ALOS.m
%             crosstrackHermiteLOS
%
% References: Fossen, T. I. and A. P. Aguiar (2024). A Uniform Semiglobal 
%             Exponential Stable Adaptive Line-of-Sight (ALOS) Guidance Law 
%             for 3-D Path Following. Automatica 163, 111556.
%             https://doi.org/10.1016/j.automatica.2024.111556
%
%             Fossen, T. I. (2023). An Adaptive Line-of-sight (ALOS) 
%             Guidance Law for Path Following of Aircraft and Marine Craft. 
%             IEEE Transactions on Control Systems Techn. 31(6), 2887-2894.
%             https://doi.org/10.1109/TCST.2023.3259819
%
% See also:   Simulink example: demoOtterUSVPathFollowingHeadingControl.slx
%
% Author:     Thor I. Fossen
% Date:       2021-04-25
% Revisions:  2023-10-14 - Added a heading autopilot and reference model
%             2024-04-01 - Added ALOS/ILOS path-following control algorithms 
%                          for straight-line paths and Hermite splines

clear ALOSpsi ILOSpsi crosstrackHermiteLOS % clear persistent variables
clearvars;

%% USER INPUTS
h  = 0.05;        % sampling time [s]
N  = 20000;		  % number of samples

% Control system flag                  
ControlFlag = 3;  % 0: PID heading autopilot, no path following
                  % 1: ALOS path-following control using straight lines 
                  %    and waypoint switching
                  % 2: ILOS path-following control using straight lines 
                  %    and waypoint switching
                  % 3: ALOS path-following control using Hermite splines

% Load condition
mp = 25;                        % payload mass (kg), max value 45 kg
rp = [0.05 0 -0.35]';           % location of payload (m)

% Ocean current
V_c = 0.3;                      % ocean current speed (m/s)
beta_c = deg2rad(30);           % ocean current direction (rad)

% Waypoints
wpt.pos.x = [0 0   150 150 -100 -100 200]';
wpt.pos.y = [0 200 200 -50  -50  250 250]';
wayPoints = [wpt.pos.x wpt.pos.y];

% LOS parameters
Delta_h = 10;                   % look-ahead distance
gamma_h = 0.001;                % ALOS adaptive gain
kappa = 0.001;                  % ILOS integral gain

% Additional parameter for straigh-line path following
R_switch = 5;                   % radius of switching circle
K_f = 0.5;                      % yaw rate observer gain

% PID heading autopilot (Nomoto gains)
T = 1;
m = 41.4;                       % m = T/K
K = T / m;

wn = 1.5;                       % closed-loop natural frequency (rad/s)
zeta = 1.0;                     % closed loop relative damping factor (-)

Kp = m * wn^2;                  % PID gains
Kd = m * (2 * zeta * wn - 1/T);
Td = Kd / Kp; 
Ti = 10 / wn;

% Reference model parameters
wn_d = 1.0;                     % natural frequency (rad/s)
zeta_d = 1.0;                   % relative damping factor (-)
r_max = deg2rad(10.0);          % maximum turning rate (rad/s)

% Otter USV input matrix
y_prop = 0.395;                 % distance from centerline to propeller (m)
k_pos = 0.0111;                 % positive Bollard, one propeller 
B = k_pos * [ 1 1               % input matrix
             y_prop  -y_prop  ];
Binv = inv(B);

% Propeller dynamics
T_n = 0.1;                      % propeller time constant (s)      
n = [0 0]';                     % intital speed n = [ n_left n_right ]'

% Inital heading, vehicle points towards next waypoint
psi0 = atan2(wpt.pos.y(2)-wpt.pos.y(1),wpt.pos.x(2)-wpt.pos.x(1));

% Initial states
eta = [0 0 0 0 0 psi0]';        % eta = [x y z phi theta psi]' 
nu  = [0 0 0 0 0 0]';           % nu  = [u v w p q r]'	 
z_psi = 0;                      % integral state
psi_d = eta(6);                 % reference model states
r_d = 0;
a_d = 0;

%% Display
disp('--------------------------------------------------------------------');
disp('MSS toolbox: Otter USV (Length = 2.0 m, Beam = 1.08 m)')  
if (ControlFlag == 0)
    disp('PID heading autopilot with reference feeforward')
elseif ControlFlag == 1
    disp(['ALOS path-following control using straight lines and ' ...
        'waypoint switching']) 
    disp(['Cirlce of acceptance: R = ',num2str(R_switch), ' m'])
    disp(['Look-ahead distance: Delta_h = ',num2str(Delta_h), ' m'])
elseif ControlFlag == 2
    disp(['ILOS path-following control using straight lines and ' ...
        'waypoint switching'])  
    disp(['Cirlce of acceptance: R =  ',num2str(R_switch), ' m'])
    disp(['Look-ahead distance: Delta_h = ',num2str(Delta_h), ' m'])
else
    disp('ALOS path-following control using Hermite splines') 
    disp(['Look-ahead distance: Delta_h = ',num2str(Delta_h), ' m'])
end
disp('--------------------------------------------------------------------');

%% MAIN LOOP
simdata = zeros(N+1,15);                   % table for simulation data

for i=1:N+1

    t = (i-1) * h;                          % time (s)

    % Measurements
    x = eta(1) + 0.01 * randn;
    y = eta(2) + 0.01 * randn;
    r = nu(6) + 0.001 * randn;
    psi = eta(6) + 0.001 * randn;

    % Guidance and control system
    if ControlFlag == 0 % heading autopilot with reference model

        % Reference model, step input
        psi_ref = psi0;
        if t > 100; psi_ref = deg2rad(0);end
        if t > 500; psi_ref = deg2rad(-90);end

        % Reference model propagation
        [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, r_max,...
            zeta_d, wn_d, h, 1);

    elseif ControlFlag == 1 % ALOS heading autopilot straight-line path following

        psi_ref = ALOSpsi(x,y,Delta_h,gamma_h,h,R_switch,wpt);
        [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

    elseif ControlFlag == 2 % ILOS heading autopilot straight-line path following

        psi_ref = ILOSpsi(x,y,Delta_h,kappa,h,R_switch,wpt);
        [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

    else % ALOS heading autopilot, Hermite spline interpolation

        undulation_factor = 0.5; % larger than 0, reduce to avoid undulation

        [psi_ref, y_e, pi_h, closestPoint, closestTangent] = ...
            crosstrackHermiteLOS(wayPoints, [x y], h,...
            undulation_factor, Delta_h, gamma_h);
        [psi_d, r_d] = LOSobserver(psi_d, r_d, psi_ref, h, K_f);

    end

   % PID heading (yaw moment) autopilot and forward thrust
   tau_X = 100;
   tau_N = (T/K) * a_d + (1/K) * r_d -...
       Kp * (ssa( psi-psi_d) +...
       Td * (r - r_d) + (1/Ti) * z_psi );
   u = Binv * [tau_X tau_N]';                   % control allocation
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

end

%% PLOTS
t = simdata(:,1); 
eta = simdata(:,2:7); 
nu  = simdata(:,8:13); 
r_d = simdata(:,14); 
psi_d = simdata(:,15); 

% Figure 1
figure(1); clf; 
hold on;
plot(eta(:,2),eta(:,1),'b');

if ControlFlag == 0 % Vehicle position

    legend('Vehicle position')

elseif ControlFlag == 1 || ControlFlag == 2 % Straigh lines and the circles

    for idx = 1:length(wayPoints(:,1))-1
        if idx == 1
            plot([wayPoints(idx,2),wayPoints(idx+1,2)],[wayPoints(idx,1),...
                wayPoints(idx+1,1)], 'r--', 'DisplayName', 'Line');
        else
            plot([wayPoints(idx,2),wayPoints(idx+1,2)],[wayPoints(idx,1),...
                wayPoints(idx+1,1)], 'r--','HandleVisibility', 'off');
        end
    end
   
    theta = linspace(0, 2*pi, 100); 
    for idx = 1:length(wayPoints(:,1))
        xCircle = R_switch * cos(theta) + wayPoints(idx,1);
        yCircle = R_switch * sin(theta) + wayPoints(idx,2);
        plot(yCircle, xCircle, 'k'); 
    end

    legend('Vehicle position','Straigh-line path','Circle of acceptance')

else % ControlFlag == 3, Hermite splines

    k = linspace(0, 1, 200);
    for i = 1:size(wayPoints, 1)-1
        P = zeros(length(k), 2);
        for j = 1:length(k)
            [P(j, :), ~] = hermiteSpline(k(j), i, wayPoints,undulation_factor);
        end
        if i == 1
            plot(P(:, 2), P(:, 1), 'r--', 'LineWidth', 3);
        else
            plot(P(:, 2), P(:, 1), 'r--', 'LineWidth', 3,...
                'HandleVisibility', 'off');
        end
    end
    
    plot(wayPoints(:, 2), wayPoints(:, 1), 'ko', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 15);
    legend('Vehicle position','Hermite spline','Waypoints')

end

xlabel('East (m)', 'FontSize', 14);
ylabel('North (m)', 'FontSize', 14);
title('North-East positions', 'FontSize', 14);
axis equal; 
hold off
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','legend'),'FontSize',14)

% Figure 2
figure(2); clf;
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

% Figure 3
figure(3); clf;
subplot(511),plot(t, sqrt(nu(:,1).^2+nu(:,2).^2));
xlabel('time (s)'),title('speed (m/s)'),grid
subplot(512),plot(t,eta(:,3),'linewidt',2)
xlabel('time (s)'),title('heave position (m)'),grid
subplot(513),plot(t,rad2deg(eta(:,4)))
xlabel('time (s)'),title('roll angle (deg)'),grid
subplot(514),plot(t,rad2deg(eta(:,5)))
xlabel('time (s)'),title('pitch angle (deg)'),grid
subplot(515),plot(t,rad2deg(unwrap(eta(:,6))),t,rad2deg(unwrap(psi_d)))
legend('\psi','\psi_d'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

