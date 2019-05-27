echo on
% SIMnavalvessel User editable script for simulation of the naval vessel under feedback control
%
% Calls:       navalvessel.m and euler2.m
%
% Author:      Thor I. Fossen
% Date:        2019-05-27
% Revisions: 

echo off
disp('Simulating navalvessel.m under PD-control with psi_ref=5 (deg) ...')

t_f = 600;   % final simulation time (sec)
h   = 0.05;   % sample time (sec)

Kp = 1e6;     % controller P-gain
Td = 20;     % controller derivative time

% initial states:
x  = [6 0 0 0 0 0 ]';   % x = [u,v p r phi psi] ]'   

% initial control inputs
tau = [1e5 0 0 0]';

% --- MAIN LOOP ---
N = round(t_f/h);               % number of samples
xout = zeros(N+1,length(x)+length(tau)+2);  % memory allocation

for i=1:N+1
    time = (i-1)*h;                   % simulation time in seconds

    r   = x(4);
    psi = x(6);
    
    % control system
    psi_ref = 5*(pi/180);                % desired heading
    tau(4) = -Kp*((psi-psi_ref)+Td*r);   % PD-controller
    
    % ship model
    [xdot,U] = navalvessel(x,tau);  % ship model
   
    % store data for presentation
    xout(i,:) = [time,x',tau',U]; 
    
    % numerical integration
    x  = euler2(xdot,x,h);             % Euler integration      
end

% time-series
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);          
p     = xout(:,4) * 180/pi;   
r     = xout(:,5) * 180/pi;
phi   = xout(:,6) * 180/pi;
psi   = xout(:,7) * 180/pi;
tau4  = xout(:,11);
U     = xout(:,12);

% plots
figure(1)
plot(t,U,'r')
grid,xlabel('time (s)'),title('speed U (m/s)'),grid

figure(2)
subplot(221),plot(t,r,'r'),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
subplot(222),plot(t,phi,'r'),xlabel('time (s)'),title('roll angle \phi (deg)'),grid
subplot(223),plot(t,psi,'r'),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
subplot(224),plot(t,tau4,'r'),xlabel('time (s)'),title('yaw moment (Nm)'),grid