function SIMtanker()
% SIMtanker is compatibel with MATLAB and GNU Octave (www.octave.org).
% This script simulates the dynamics of a large tanker, length 304.8 m, 
% under feedback control. 
%
% Dependencies:
%   tanker     - Nonlinear tanker model
%   euler2     - Euler's integrations method
%
% Reference: 
%   W. B. Van Berlekom and T. A. and Goddard (1972). Maneuvering of Large
%     Tankers, Transaction of SNAME, 80:264-298
%
% Author:    Thor I. Fossen
% Date:      2024-05-02
% Revisions:
%   None

close all;
clearvars;

%% USER INPUTS
t_f = 3000;             % Final simulation time (sec)
h   = 0.1;              % Sample time (sec)

water_depth = 100;      % Water depth must be larger than draft 18.5 m
psi_ref = deg2rad(5);   % Desired heading
wn = 0.1;               % Closed-loop natural frequency (rad/s)
Kp = 10;                % Controller P gain
Td = 10;                % Controller derivative time

% initial states:
x = [7 0 0 0 0 0 0 60]';   % x = [ u v r x y psi delta n ]'
psi_d = 0;

%% MAIN LOOP
N = round(t_f/h);                       % Number of samples
simdata = zeros(N+1,length(x)+2);       % Memory allocation
                
for i=1:N+1

    time = (i-1) * h;                   % Simulation time in seconds

    r   = x(3);
    psi = x(6);
    
    % Control system (constant thrust + PD heading controller)
    if time > 1200
        psi_d = exp(-h*wn) * psi_d + (1 - exp(-h*wn)) * psi_ref;
    end
    delta_c = -Kp * ( ssa(psi-psi_d) + Td * r );  % PD controller
    n_c = 70;
    
    % Ship model
    [xdot,U] = tanker(x,[-delta_c n_c water_depth]);   % Tanker model
   
    % Store data for presentation
    simdata(i,:) = [time, x', U]; 
    
    % Euler's integration method (k+1)
    x = euler2(xdot,x,h);                         

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
r     = rad2deg(simdata(:,4));   
x     = simdata(:,5);
y     = simdata(:,6);
psi   = rad2deg(simdata(:,7));
delta = rad2deg(simdata(:,8));
n     = simdata(:,9);
U     = simdata(:,10);

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
subplot(221),plot(t,r),xlabel('Time (s)')
title('Yaw rate r (deg/s)'),grid
subplot(222),plot(t,U),xlabel('Time (s)')
title('Speed U (m/s)'),grid
subplot(223),plot(t,psi),xlabel('Time (s)')
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
    'Mass (no payload)', '55.0 kg',...
    'Volume displacement', '220 000 m3',...
    'Nominal speed', '16 knots',...
    'Max rudder angle', '10 deg',...
    'Max propeller speed', '80 RPM'};
displayVehicleData('Esso 190 000 dwt Tanker', vesselData, 'tanker.jpg', 3);

end
