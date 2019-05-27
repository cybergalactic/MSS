echo on
% SIMcontainer User editable script for simulation of the container ship under feedback control
%              Both the linear model Lcontainer.m and nonlinear model container.m are simulated.
%
% Calls:       container.m, Lcontainer.m and euler2.m
%
% Author:      Thor I. Fossen
% Date:        2018-07-21
% Revisions: 

echo off
disp('Simulating container.m and Lcontainer.m under PD-control with psi_ref=5 (deg) ...')

t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)

Kp = 1;      % controller P-gain
Td = 10;     % controller derivative time

% initial states:
x  = [7 0 0 0 0 0 0 0 0 70]';   % x = [u v r x y psi p phi delta n ]'
x2 = [7 0 0 0 0 0 0 0 0]';     

% --- MAIN LOOP ---
N = round(t_f/h);               % number of samples
xout = zeros(N+1,length(x)+2);  % memory allocation
xout2 = zeros(N+1,length(x2)+2);  % memory allocation

for i=1:N+1
    time = (i-1)*h;                   % simulation time in seconds

    r   = x(3);
    psi = x(6);
    
    % control system
    psi_ref = 5*(pi/180);                % desired heading
    delta_c = -Kp*((psi-psi_ref)+Td*r);  % PD-controller
    n_c = 70;
    
    % ship model
    [xdot,U] = container(x,[delta_c n_c]);  % ship model, see .../gnc/VesselModels/
    [xdot2,U2] = Lcontainer(x2,delta_c);
   
    % store data for presentation
    xout(i,:) = [time,x',U]; 
    xout2(i,:) = [time,x2',U2]; 
    
    % numerical integration
    x  = euler2(xdot,x,h);             % Euler integration
    x2 = euler2(xdot2,x2,h);        
end

% time-series
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);          
r     = xout(:,4)*180/pi;   
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7)*180/pi;
p     = xout(:,8)*180/pi;
phi   = xout(:,9)*180/pi;
delta = xout(:,10)*180/pi;
n     = xout(:,11);
U     = xout(:,12);

t2     = xout2(:,1);
u2     = xout2(:,2); 
v2     = xout2(:,3);          
r2     = xout2(:,4)*180/pi;   
x2     = xout2(:,5);
y2     = xout2(:,6);
psi2   = xout2(:,7)*180/pi;
p2     = xout2(:,8)*180/pi;
phi2   = xout2(:,9)*180/pi;
delta2 = xout2(:,10)*180/pi;
U2     = xout2(:,11);

% plots
figure(1)
plot(y,x,'r',y2,x2,'g')
grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position')
legend('nonlinear model','linear model')

figure(2)
subplot(221),plot(t,r,'r',t,r2,'g'),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
legend('nonlinear model','linear model')
subplot(222),plot(t,phi,'r',t,phi2,'g'),xlabel('time (s)'),title('roll angle \phi (deg/s)'),grid
legend('nonlinear model','linear model')
subplot(223),plot(t,psi,'r',t,psi2,'g'),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
legend('nonlinear model','linear model')
subplot(224),plot(t,delta,'r',t,delta2,'g'),xlabel('time (s)'),title('rudder angle \delta (deg)'),grid
legend('nonlinear model','linear model')