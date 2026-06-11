function SIMnavalvessel()
% SIMnavalvessel is compatibel with MATLAB and GNU Octave (www.octave.org). 
% This script simulates a naval vessel under PD heading control.
%
% Dependencies:  
%   navalvessel.m   - Naval vessel dynamics (Blanke and Christensen, 1993).
%
% Reference: 
%   M. Blanke and A. Christensen (1993). Rudder-roll damping autopilot 
%   robustness to sway-yaw-roll couplings. 10th Ship Control 
%   Systems Symposium, Ottawa, Canada.
%
% Author:      Thor I. Fossen
% Date:        2019-05-27
% Revisions:
%   2024-03-27 : Added animation of the ship North-East positions.
%   2024-04-22 : Enhanced compatibility with GNU Octave.

clear animateShip   % clear the persistent animation variables
clearvars;

T_final = 600;	    % Final simulation time (s)
h   = 0.05;         % Sample time (s)

% PD controller
Kp = 1e6;           % Controller proportional gain
Td = 1;            % Controller derivative time

% Reference signal
w_n = 0.01;             % Low-pass filter cut-off frequency (rad/s)
psi_ref = deg2rad(20);  % Desired heading

% Initial states
x  = [6 0 0 0 0 0 ]';   % x = [u v p r phi psi] ]'   
eta = [0 0 0]';         % initial position expressed in NED
xf_psi_d = eta(3);

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% Display simulation options
displayControlMethod();

%% MAIN LOOP
simdata = zeros(nTimeSteps,11);       % Preallocate table 

for i=1:nTimeSteps

    % Measurements
    r   = x(4) + 0.001 * randn;
    psi = x(6) + 0.001 * randn;

    % Desired heading
    [xf_psi_d, psi_d] = lowPassFilter(xf_psi_d, psi_ref, w_n, h);  
    
    % Control system
    tauX = 1e5;                                  % Thrust
    tauN = -Kp * ( ssa(psi - psi_d) + Td * r );  % PD heading controller
    
    % Ship dynamics
    tau = [tauX 0 0 tauN]';
    [xdot,U] = navalvessel(x,tau);       
   
    % Store data for presentation
    simdata(i,:) = [x', tauN, U, eta(1:2)', psi_d]; 
    
    % Numerical integration
    x = euler2(xdot,x,h);                 
    eta = eta + h * Rzyx(0,0,psi) * [x(1) x(2) x(4)]'; 

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

u     = simdata(:,1); 
v     = simdata(:,2);          
p     = rad2deg(simdata(:,3));   
r     = rad2deg(simdata(:,4));
phi   = rad2deg(simdata(:,5));
psi   = rad2deg(simdata(:,6));
tauN  = simdata(:,7);
U     = simdata(:,8);
x     = simdata(:,9);
y     = simdata(:,10);
psi_d = rad2deg(simdata(:,11));

% Plot and animation of the North-East positions
figure(1)
shipSize = 0.2;
set(gcf, 'Position', [1, 1, 0.4*scrSz(3), scrSz(4)]);
animateShip(x,y,shipSize,'b-',1);

figure(2)
plot(t,U,'b')
xlabel('time (s)'),title('Speed U (m/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(3)
subplot(221)
plot(t,r,'b'),xlabel('time (s)'),title('Yaw rate r (deg/s)'),grid
subplot(222)
plot(t,phi,'b'),xlabel('time (s)'),title('Roll angle \phi (deg)'),grid
subplot(223)
plot(t,psi,'b',t, psi_d,'r'),xlabel('time (s)'),
title('Yaw angle \psi (deg)'),grid
subplot(224)
plot(t,tauN,'b'),xlabel('time (s)'),title('Yaw control moment (Nm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', '51.5 m', ...
    'Beam', '8.6 m', ...
    'Draft', '2.3 m', ...
    'Mass', '362 tonnes', ...
    'Cruise speed', '8.23 m/s'};
displayVehicleData('Multipurpose Naval Vessel', vesselData, 'nvessel.jpg', 4);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Multipurpose Naval Vessel');
    disp('Norrbin (1963) nonlinear model');    
    disp('Heading autopilot: PD control law')
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end

