echo on
% SIMDEMO1    User editable script for simulation of the 
%             Mariner class vessel under feedback control
%
% Calls:      mariner.m
%             euler2.m
%
% Author:     Thor I. Fossen
% Date:       19 Ju 2001
% Revisions: 

echo off 
disp('Simulating mariner.m under PD-control with psi_ref=5 (deg) ...')

t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)

Kp = 1;      % controller P-gain
Td = 10;     % controller derivative time

% initial states:  x = [ u v r x y psi delta ]' 
x = zeros(7,1);   

%% MAIN LOOP
N = round(t_f/h);               % number of samples
xout = zeros(N+1,length(x)+2);  % memory allocation

for i=1:N+1
    time = (i-1)*h;                   % simulation time in seconds

    r   = x(3);
    psi = x(6);
    
    % control system
    psi_ref = 5*(pi/180);              % desired heading
    delta = -Kp*((psi-psi_ref)+Td*r);  % PD-controller

    % ship model
    [xdot,U] = mariner(x,delta);       % ship model, see .../gnc/VesselModels/
    
    % store data for presentation
    xout(i,:) = [time,x',U]; 
    
    % numerical integration
    x = euler2(xdot,x,h);             % Euler integration
end

%% PLOTS
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);          
r     = xout(:,4)*180/pi;   
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7)*180/pi;
delta = xout(:,8)*180/pi;
U     = xout(:,9);

figure(1)
plot(y,x),grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(221),plot(t,r),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
subplot(222),plot(t,U),xlabel('time (s)'),title('speed U (m/s)'),grid
subplot(223),plot(t,psi),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
subplot(224),plot(t,delta),xlabel('time (s)'),title('rudder angle \delta (deg)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
