function zigzag6dof(vehicle,x,ui,t_final,t_rudderexecute,h,maneuver)
% zigzag6dof(vehicle,x,ui,t_final,t_rudderexecute,h,maneuver) performs the 
% zigzag maneuver for 6-DOF models, see ExZigZag.m
%
% Inputs:
% 'vehicle':       vehicle model. Compatible with the Remmus 100 AUV located
%                  under .../VESSELS
% x:               initial state vector for vehicle model
% ui:              [delta,:] where delta=0 and the other values are nonzero if any
% t_final:         final simulation time
% t_rudderexecute: time for which control input is activated
% h:               sampling time
% maneuver         [rudder angle, heading angle]. Default 20-20 deg, i.e. 
%                  maneuver = [20, 20] rudder is changed to maneuver(1) when 
%                  the heading angle is larger than maneuver(2)
%
% Author:          Thor I. Fossen
% Date:            18 Nov 2023
% Revisions: 

if nargin>7 || nargin < 6, error('number of inputs must be 6 or 7'); end
if t_final < t_rudderexecute
    error('t_final must be larger than t_rudderexecute'); 
end
if nargin == 6, maneuver = [20,20]; end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,15);               % memory allocation

disp('Simulating...')

u_vehicle = ui;

for i=1:N+1
    time = (i-1)*h;
    
    psi = rad2deg( x(12) );
    r = x(6);
    
    if round(time) == t_rudderexecute 
        u_vehicle(1) = deg2rad( maneuver(1) ); 
    end
    
    if round(time) > t_rudderexecute
        if ( psi >= maneuver(2) && r > 0 )
            u_vehicle(1) = -deg2rad( maneuver(1) ); 
        elseif (psi <= -maneuver(2) && r < 0 )
            u_vehicle(1) = deg2rad( maneuver(1) );            
        end   
    end
 
    [xdot,U] = feval(vehicle,x,u_vehicle);          % vehicle model
    
    xout(i,:) = [time,x(1:12)',U,u_vehicle(1)];  
    
    x = euler2(xdot,x,h);                     % Euler's method
end

% time-series
t       = xout(:,1);
u       = xout(:,2); 
v       = xout(:,3); 
w       = xout(:,4); 
p       = rad2deg( xout(:,5) );   
q       = rad2deg( xout(:,6) );   
r       = rad2deg( xout(:,7) );   
x       = xout(:,8);
y       = xout(:,9);
z       = xout(:,10);
phi     = rad2deg( xout(:,11) );
theta   = rad2deg( xout(:,12) );
psi     = rad2deg( xout(:,13) );
U       = xout(:,14);
delta_c = rad2deg( xout(:,15) );

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
