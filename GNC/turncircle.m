function [t,u,v,r,x,y,psi,U] = turncircle(ship,x,ui,t_final,t_rudderexecute,h)
% TURNCIRCLE  [t,u,v,r,x,y,psi,U] = turncircle(ship,x,ui,t_final,t_rudderexecute,h)
% computes the turning circle maneuvering indexes, see ExTurnCircle.m
%
% Inputs:
%  'ship':          Ship model, compatible with the models under /MSS/Vessels/.
%  x:               Initial state vector for ship model
%  ui:              ui = [delta,:] where delta is the rudder command at 
%                   time = t_rudderexecute
%  t_final:         Final simulation time
%  t_rudderexecute: Time control input is activated
%  h:               Sampling time
%
% Outputs:
%  t:               Time vector
%  u,v,r,x,y,psi,U: Time series
%
% Author:    Thor I. Fossen
% Date:      18 Jul 2001
% Revisions: 25 Nov 2002 - Expression for Nrudder was corrected, included
%                          plots for rudder execute, 90 deg heading angle

if nargin~=6, error('number of inputs must be 6'); end
if t_final<t_rudderexecute, error('t_final must be larger than t_rudderexecute'); end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,8);                % memory allocation
store1 = 1; store2 = 1;             % logical variables (0,1)

disp('Simulating...')

for i=1:N+1

    time = (i-1)*h;
    
    if round(abs(x(6))*180/pi)>=90 & store1==1
        transfer=x(5);   % transfer at 90 deg
        advance =x(4);   % advance at 90 deg
        store1 = 0;
    end
    
    if round(abs(x(6))*180/pi)>=180 & store2==1
        tactical=x(5);   % tactical diameter at 180 deg
        store2 = 0;
    end
    
    u_ship = ui;
    if round(time) < t_rudderexecute
       u_ship(1) = 0;   % zero rudder angle
    end     
    
    [xdot,U] = feval(ship,x,u_ship);       % ship model
    
    xout(i,:) = [time,x(1:6)',U];  
    
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

Nrudder = round(t_rudderexecute/h);

% turning radius, tactical diameter, advance and transfer
disp(' ')
fprintf('Rudder execute (x-coordinate)          : %4.0f m\n',abs(x(Nrudder)))
fprintf('Steady turning radius                  : %4.0f m\n',U(N+1)/abs(r(N+1)*pi/180))
fprintf('Maximum transfer                       : %4.0f m\n',abs(max(abs(y))))
fprintf('Maximum advance                        : %4.0f m\n',abs(max(abs(x))-x(Nrudder)))      
fprintf('Transfer at 90 (deg) heading           : %4.0f m\n',abs(transfer))    
fprintf('Advance at 90 (deg) heading            : %4.0f m\n',abs(advance-x(Nrudder)))        
fprintf('Tactical diameter at 180 (deg) heading : %4.0f m\n',abs(tactical))

% plots
figure(1)
plot(x,y,x(Nrudder),y(Nrudder),'linewidth',2), hold on
plot(x(Nrudder),y(Nrudder),'*r',advance,transfer,'or'), hold off
grid,axis('equal'),xlabel('x-position'),ylabel('y-position')
title('Turning circle (* = rudder execute, o = 90 deg heading)')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(211),plot(t,r),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
subplot(212),plot(t,U),xlabel('time (s)'),title('speed U (m/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)


