function [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute,h,maneuver)
% ZIGZAG      [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute,h,maneuver)
%             performs the zig-zag maneuver, see ExZigZag.m
%
% Inputs :
% 'ship'          = ship model. Compatible with the models under .../gnc/VesselModels/
% x               = initial state vector for ship model
% ui              = [delta,:] where delta=0 and the other values are non-zero if any
% t_final         = final simulation time
% t_rudderexecute = time control input is activated
% h               = sampling time
% maneuver        = [rudder angle, heading angle]. Default 20-20 deg that is: maneuver = [20, 20] 
%                    rudder is changed to maneuver(1) when heading angle is larger than maneuver(2)
%
% Outputs :
% t               = time vector
% u,v,r,x,y,psi,U = time series
%
% Author:    Thor I. Fossen
% Date:      22th July 2001
% Revisions: 15th July 2002, switching logic has been modified to handle arbitrarily maneuvers

if nargin>7 | nargin<6, error('number of inputs must be 6 or 7'); end
if t_final<t_rudderexecute, error('t_final must be larger than t_rudderexecute'); end
if nargin==6, maneuver = [20,20]; end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,9);                % memory allocation

disp('Simulating...')

u_ship=ui;

for i=1:N+1,
    time = (i-1)*h;
    
    psi = x(6)*180/pi;
    r   = x(3);
    
    if round(time)==t_rudderexecute, 
        u_ship(1)=maneuver(1)*pi/180; 
    end
    
    if round(time) > t_rudderexecute, 
        if (psi>=maneuver(2) & r>0),
            u_ship(1) = -maneuver(1)*pi/180; 
        elseif (psi<=-maneuver(2) & r<0),
            u_ship(1) = maneuver(1)*pi/180;            
        end   
    end
 
    [xdot,U] = feval(ship,x,u_ship);       % ship model
    
    xout(i,:) = [time,x(1:6)',U,u_ship(1)];  
    
    x = euler2(xdot,x,h);                     % Euler integration
end

% time-series
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);         
r     = xout(:,4)*180/pi; 
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7)*180/pi;
U     = xout(:,8);
delta_c = xout(:,9)*180/pi;

% plots
figure(1)
plot(x,y,'linewidth',2),grid,axis('equal'),xlabel('x-position'),ylabel('y-position')
title('Zig-zag test')
figure(2)
subplot(211),plot(t,psi,'linewidth',2)
hold on
plot(t,delta_c,'r')
hold off
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
legend('\psi','\delta_c')
subplot(212),plot(t,U),xlabel('time (s)'),title('speed U (m/s)'),grid
