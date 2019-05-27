echo on
% SIMmariner  User editable script for simulation of the 
%             mariner class vessel under feedback control
%
% Calls:      mariner.m and euler2.m
%
% Author:     Thor I. Fossen
% Date:       2018-07-21
% Revisions: 

echo off 
disp('Simulating mariner.m under PD-control with psi_ref=5 (deg) ...')

t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)

Kp = 1;      % controller P-gain
Td = 10;     % controller derivative time

% initial states:  x = [ u v r x y psi delta ]' 
x = zeros(7,1);   

% --- MAIN LOOP ---
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
    [xdot,U] = mariner(x,delta);       % ship model
    
    % store data for presentation
    xout(i,:) = [time,x',U]; 
    
    % numerical integration
    x = euler2(xdot,x,h);             % Euler integration
end

% time-series
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);          
r     = xout(:,4)*180/pi;   
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7)*180/pi;
delta = xout(:,8)*180/pi;
U     = xout(:,9);

% plots
figure(1)
plot(y,x),grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position')

figure(2)
subplot(221),plot(t,r),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
subplot(222),plot(t,U),xlabel('time (s)'),title('speed U (m/s)'),grid
subplot(223),plot(t,psi),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
subplot(224),plot(t,delta),xlabel('time (s)'),title('rudder angle \delta (deg)'),grid
