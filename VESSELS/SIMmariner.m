function SIMmariner()
% SIMmariner is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script simulates the dynamic behavior of a Mariner-Class Cargo Vessel, 
% length 160.93 m, under PID heading control, and waypoint path-following 
% control using a course autopilot (Fossen 2022). The Speed Over Ground 
% (SOG) and Course over Ground (COG) are estimated during path following 
% using the 5-state Extended Kalman Filter (EKF) by Fossen and Fossen (2021). 
%
% Dependencies:
%   mariner.m           - Vessel dynamics.  
%   EKF_5states.m       - EKF for estimation of SOG and COG.
%   refModel.m          - Reference model for autopilot systems.
%   LOSchipsi.m         - LOS guidance algorithm for path following.
%   LOSobserver.m       - Observer for LOS guidance. 
%   controlMethods.m    - Menu for choosing control law.
%
% Simulink Models:
%   demoMarinerPathFollowingCourseControl.slx
%
% References: 
%   T. I. Fossen (2022). Line-of-sight Path-following Control utilizing an 
%      Extended Kalman Filter for Estimation of Speed and Course over Ground 
%      from GNSS Positions. Journal of Marine Science and Technology 27, 
%      pp. 806–813.
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%       Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%   S. Fossen and T. I. Fossen (2021). Five-state Extended Kalman 
%      Filter for Estimation of Speed Over Ground (SOG), Course Over Ground 
%      (COG) and Course Rate of Unmanned Surface Vehicles (USVs): 
%      Experimental Results. Sensors 21(23). 
%
% Author:     Thor I. Fossen
% Date:       2018-07-21
% Revisions:
%   2024-03-27 : Added animation of the ship's North-East positions.
%   2024-04-19 : Enhanced compatibility with GNU Octave.
%   2024-05-16 : Added an option for LOS path-following course controll.

clearvars; 
close all;
clear LOSchi EKF_5states  % Clear persistent variables in functions

t_f = 3000;               % Final simulation time (sec)
h  = 0.05;                % Sampling time [s]
Z = 2;                    % GNSS measurement frequency (2 times slower)

% Waypoints
wpt.pos.x = [0 2000 5000 3000 6000 10000]';
wpt.pos.y = [0 0 5000  8000 12000 12000]';
wayPoints = [wpt.pos.x wpt.pos.y];

% LOS parameters
Delta_h = 500;                   % Look-ahead distance
R_switch = 400;                  % Radius of switching circle
K_f = 0.2;                       % LOS observer gain

% Initial heading, vehicle points towards next waypoint
psi0 = atan2(wpt.pos.y(2) - wpt.pos.y(1), wpt.pos.x(2) - wpt.pos.x(1));

% PID pole placement algorithm (Fossen 2021, Section 15.3.4)
wn = 0.05;                       % Closed-loop natural frequency
T = 107.3;                       % Nomoto time constant
K = 0.185;                       % Nomoto gain constant
Kp = (T/K) * wn^2;               % Proportional gain
Td = T/(K*Kp) * (2*1*wn - 1/T);  % Derivative time constant
Ti = 10 / wn;                    % Integral time constant

% Reference model specifying the heading autopilot closed-loop dynamics
wn_d = 0.1;                      % Natural frequency (rad/s)
zeta_d = 1.0;                    % Relative damping factor (-)
r_max = deg2rad(1.0);            % Maximum turning rate (rad/s)

% Initial states
x_hat = zeros(5,1);              % xhat = [ x, y, U, chi, oemga_chi ]'
x = zeros(7,1);                  % x = [ u v r x y psi delta ]'
U0 = 7.7175;                     % Nominal speed
e_int = 0;                       % Autopilot integral state 
delta_c = 0;                     % Initial rudder angle command
psi_d = 0;                       % Initial desired heading angle
chi_d = 0;                       % Initial desired course angle
omega_chi_d = 0;                 % Initial desired course rate
r_d = 0;                         % Initial desired rate of turn
a_d = 0;                         % Initial desired acceleration

% Choose control method and display simulation options
methods = {'PID heading autopilot, no path following',...
           'LOS path-following using a course autopilot with waypoint switching'};
ControlFlag = controlMethod(methods);
displayControlMethod(ControlFlag, R_switch, Delta_h);

%% MAIN LOOP
N = round(t_f/h);                    % Number of samples
simdata = zeros(N+1,17);             % Memory allocation

for i=1:N+1

    t = (i-1) * h;                   % Simulation time in seconds

    % Measurements with measurement noise    
    r    = x(3) + 0.0001 * randn;
    xpos = x(4) + 0.01 * randn;
    ypos = x(5) + 0.01 * randn;
    psi  = x(6) + 0.0001 * randn;

    % EKF estimates used for path-following control
    U_hat = x_hat(3);
    chi_hat = x_hat(4);
    omega_chi_hat = x_hat(5);

    % Guidance and control system
    switch ControlFlag

        case 1  
            % PID heading autopilot        
            psi_ref = psi0; % Reference model, step input adjustments
            if t > 100;  psi_ref = deg2rad(30); end
            if t > 1000; psi_ref = deg2rad(-30); end

            % PID heading autopilot
            e = ssa(psi - psi_d);
            delta_PID = (T/K) * a_d + (1/K) * r_d ...         % Feedforward
               -Kp * ( e + Td * (r - r_d) + (1/Ti) * e_int ); % PID                             

            % Reference model propagation
            [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, psi_ref, ...
                r_max, zeta_d, wn_d, h, 1);

        case 2  
            % LOS course autopilot for straight-line path following
            chi_ref = LOSchi(xpos, ypos, Delta_h, R_switch, wpt);

            % LOS observer for estimation of chi_d and omega_chi_d
            [chi_d, omega_chi_d] = LOSobserver(...
                chi_d, omega_chi_d, chi_ref, h, K_f);

            omega_chi_d = sat(omega_chi_d, deg2rad(1)); % Max value

            % PID course autopilot
            e = ssa(chi_hat - chi_d);
            delta_PID = (1/K) * omega_chi_d ...               % Feedforward
               -Kp * ( e + Td * (omega_chi_hat - omega_chi_d) ... % PID
               + (1/Ti) * e_int );                 

    end
    
    delta_c = sat(delta_c, deg2rad(40));             % Maximum rudder angle

    % Ship dynamics
    [xdot,U] = mariner(x,delta_c);     
    
    % Store data for presentation
    simdata(i,:) = [t, x', U, psi_d, r_d, chi_d, omega_chi_d, delta_c, ...
        U_hat, chi_hat, omega_chi_hat]; 
    
    % Numerical integration
    x = x + h * xdot;                             % Euler's method
    e_int = e_int + h * e;
    delta_c = delta_c + h * (delta_PID - delta_c) / 1.0;

    % Propagation of the EKF states
    x_hat = EKF_5states(xpos, ypos, h, Z, 'NED', ...
        100*diag([0.1,0.1]), 1000*diag([1 1]), 0.00001, 0.00001);

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

% Simdata(i,:) = [t, x', U, psi_d, r_d, chi_d, omega_chi_d, delta_c]
t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));
delta = -rad2deg(simdata(:,8));     % delta = -delta_r (physical angle)
U     = simdata(:,9);
psi_d = rad2deg(simdata(:,10));
r_d   = rad2deg(simdata(:,11));
chi_d = rad2deg(simdata(:,12));
omega_chi_d = rad2deg(simdata(:,13));
delta_c = rad2deg(simdata(:,14));
U_hat = simdata(:,15);
chi_hat = rad2deg(simdata(:,16));
omega_chi_hat = rad2deg(simdata(:,17));

chi = psi + atan2(v, U0 + u);

% Plot and animation of the North-East positions
figure(1)
set(gcf,'Position',[1,1, 1.0*scrSz(4),1.0*scrSz(4)],'Visible','off');
hold on;
plot(y,x,'b');  % Plot vehicle position

% Control method specific plots
if ControlFlag == 1  % Heading autopilot
    legend('Vessel position');
else   % Path-following controller, straight line and circles
    plotStraightLinesAndCircles(wayPoints, R_switch);
end

xlabel('East (m)', 'FontSize', 14);
ylabel('North (m)', 'FontSize', 14);
title('North-East Positions (m)', 'FontSize', 14);
axis equal;
grid on;
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','legend'),'FontSize',12);
set(1,'Visible', 'on');  % Show figure

figure(2)
subplot(221)
if ControlFlag == 1  
    % Heading autopilot
    plot(t,r_d,t,r)
    xlabel('Time (s)')
    title('Yaw rate (deg/s)')
    legend('Desired','True')
else
    % Course autopilot
    plot(t,omega_chi_hat,t,omega_chi_d)
    xlabel('Time (s)')
    title('Course rate (deg/s)')
    legend('Estimated','Desired')
end
grid
subplot(222)
plot(t,U_hat,t,U)
xlabel('Time (s)')
title('Speed (m/s)')
legend('Estimated','True')
grid
subplot(223)
if ControlFlag == 1  
    % Heading autopilot
    plot(t,psi_d,t,psi)
    xlabel('Time (s)')
    title('Yaw angle (deg)')
    legend('Desired','True')
else
    % Course autopilot
    plot(t,chi_hat,t,chi,t,chi_d,'k')
    xlabel('Time (s)')
    title('Course angle (deg)')
    legend('Estimated','True','Desired')
end
grid
subplot(224)
plot(t,delta,t,delta_c)
xlabel('Time (s)')
title('Rudder angle (deg)')
legend('Actual','Commanded')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', '160.93 m', ...
    'Mass', '17 045 tonnes', ...
    'Max speed', '7.71 m/s', ...
    'Max rudder angle', '40 deg'};
displayVehicleData('Mariner-Class Cargo Vessel', vesselData, 'mariner.jpg', 3);

end


function plotStraightLinesAndCircles(wayPoints, R_switch)
    legendLocation = 'northeast';
    
    % Plot straight lines and circles for straight-line path following
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

    legend('Vessel position','Straight-line path','Circle of acceptance',...
        'Location',legendLocation);
end


%% DISPLAY CONTROL METHOD
function displayControlMethod(ControlFlag, R_switch, Delta_h)
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Mariner-Class Cargo Vessel');
    disp('Five-state EKF for estimation of SOG and COG');
    switch ControlFlag
        case 1
            disp('PID course autopilot with reference feedforward');
        case 2
            disp(['LOS path-following using straight lines and ' ...
                'waypoint switching']);
            disp(['Circle of acceptance: R = ', num2str(R_switch), ' m']);
            disp(['Look-ahead distance: Delta_h = ', num2str(Delta_h), ' m']);
    end
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end
