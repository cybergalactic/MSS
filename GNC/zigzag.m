function [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute, ...
    h,maneuver)
% [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute,h,maneuver)
% performs the zigzag maneuver, see ExZigZag.m
%
% Inputs :
% 'ship':          ship model. Compatible with the models under .../VESSELS
% x:               initial state vector for ship model
% ui:              [delta,:] where delta=0 and the other values are nonzero if any
% t_final:         final simulation time
% t_rudderexecute: time for which control input is activated
% h:               sampling time
% maneuve          [rudder angle, heading angle]. Default 20-20 deg, i.e. 
%                  maneuver = [20, 20] rudder is changed to maneuver(1) when 
%                  heading angle is larger than maneuver(2)
%
% Outputs :
% t               = time vector
% u,v,r,x,y,psi,U = time series
%
% Author:    Thor I. Fossen
% Date:      22th July 2001
% Revisions: 15th July 2002, switching logic has been modified to handle 
%                            arbitrarily maneuvers

if nargin>7 || nargin < 6, error('number of inputs must be 6 or 7'); end
if t_final < t_rudderexecute
    error('t_final must be larger than t_rudderexecute'); 
end
if nargin == 6, maneuver = [20,20]; end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,9);                % memory allocation

disp('Simulating...')

u_ship=ui;

for i=1:N+1
    time = (i-1)*h;
    
    psi = rad2deg( x(6) );
    r   = x(3);
    
    if round(time) == t_rudderexecute 
        u_ship(1) = deg2rad( maneuver(1) ); 
    end
    
    if round(time) > t_rudderexecute
        if ( psi >= maneuver(2) && r > 0 )
            u_ship(1) = -deg2rad( maneuver(1) ); 
        elseif (psi <= -maneuver(2) && r < 0 )
            u_ship(1) = deg2rad( maneuver(1) );            
        end   
    end
 
    [xdot,U] = feval(ship,x,u_ship);          % ship model
    
    xout(i,:) = [time,x(1:6)',U,u_ship(1)];  
    
    x = euler2(xdot,x,h);                     % Euler's method
end

% time-series
t       = xout(:,1);
u       = xout(:,2); 
v       = xout(:,3);         
r       = rad2deg( xout(:,4) ); 
x       = xout(:,5);
y       = xout(:,6);
psi     = rad2deg( xout(:,7) );
U       = xout(:,8);
delta_c = rad2deg( xout(:,9) );

% plots
figure(1)
plot(y,x),grid,axis('equal'),xlabel('East'),ylabel('North')
title('Zigzag test')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(211),plot(t,psi,t,delta_c,'r')
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
legend('\psi','\delta_c')
subplot(212),plot(t,U),xlabel('time (s)'),title('speed U (m/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
