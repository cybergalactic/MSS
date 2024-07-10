function SIMtanker()
% SIMtanker is compatibel with MATLAB and GNU Octave (www.octave.org).
% This script simulates the dynamics of a large tanker, length 304.8 m, 
% under feedback control. 
%
% Dependencies:
%   tanker        - Nonlinear tanker model
%   lowPassdilter - Low-pass filter
%
% Reference: 
%   W. B. Van Berlekom and T. A. and Goddard (1972). Maneuvering of Large
%     Tankers, Transaction of SNAME, 80:264-298
%
% Author:    Thor I. Fossen
% Date:      2024-05-02
% Revisions:

close all;
clearvars;

%% USER INPUTS
% Define simulation parameters
T_final = 2000;	                % Final simulation time (s)
h = 0.1;                        % Sampling time (s)

water_depth = 100;      % Water depth must be larger than draft 18.5 m
psi_ref = deg2rad(5);   % Desired heading
wn = 0.1;               % Closed-loop natural frequency (rad/s)
Kp = 10;                % Controller P gain
Td = 10;                % Controller derivative time

% Initial states:
x = [7 0 0 0 0 0 0 60]';   % x = [ u v r x y psi delta n ]'
psi_d = 0;

% Display simulation options
displayControlMethod();

%% MAIN LOOP
t = 0:h:T_final;                         % Time vector
simdata = zeros(length(t), length(x)+2); % Pre-allocate matrix for efficiency
                
for i=1:length(t)

    r   = x(3);
    psi = x(6);
    
    % Control system (constant thrust + PD heading controller)
    if t(i) > 600
        psi_d = lowPassFilter(psi_d, psi_ref, wn, h);
    end
    delta_c = -Kp * ( ssa(psi-psi_d) + Td * r );  % PD controller
    n_c = 70;
    
    % Ship model
    [xdot,U] = tanker(x,[-delta_c n_c water_depth]);   % Tanker model
   
    % Store data for presentation
    simdata(i,:) = [x', U, psi_d]; 
    
    % Euler's integration method (k+1)
    x = euler2(xdot,x,h);                         

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

u     = simdata(:,1); 
v     = simdata(:,2);          
r     = rad2deg(simdata(:,3));   
x     = simdata(:,4);
y     = simdata(:,5);
psi   = rad2deg(simdata(:,6));
delta = rad2deg(simdata(:,7));
n     = simdata(:,8);
U     = simdata(:,9);
psi_d = rad2deg(simdata(:,10));

% North-East positions
figure(1); 
set(gcf, 'Position', [1, 1, 0.5*scrSz(3), scrSz(4)]);
plot(y,x)
grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position (m)')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

% Ship speed, yaw rate, yaw angle, rudder angle, and propeller speed
figure(2); 
if ~isoctave; set(gcf,'Position',[scrSz(3)/2,1,scrSz(3)/2.5,scrSz(4)]); end
subplot(221)
plot(t,r),xlabel('Time (s)')
title('Yaw rate r (deg/s)'),grid
subplot(222),plot(t,U),xlabel('Time (s)')
title('Speed U (m/s)'),grid
subplot(223)
plot(t,psi,t,psi_d),xlabel('Time (s)')
title('Yaw angle \psi (deg)'),grid
subplot(224)
plot(t,delta),xlabel('Time (s)')
hold on
plot(t,n),xlabel('Time (s)')
hold off
legend('Rudder angle (deg)','Propeller speed (RPM)')
title('Control inputs'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', '304.8 m',...
    'Beam', '47.17 m',...
    'Draft to design waterline ', '18.46 m',...    
    'Mass', '225 500 tonnes',...
    'Volume displacement', '220 000 m3',...
    'Nominal speed', '16 knots',...
    'Max rudder angle', '10 deg',...
    'Max propeller speed', '80 RPM'};
displayVehicleData('Esso 190 000 dwt Tanker', vesselData, 'tanker.jpg', 3);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Esso 190 000 dwt Tanker');
    disp('PD heading autopilot');
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end

