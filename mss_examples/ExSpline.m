% ExSpline   Cubic Hermite and spline interpolation of way-points 
% Author:    Thor I. Fossen
% Date:      6 July 2002
% Revisions: 

% way-point database
wpt.pos.x   = [0 100 500 700 1000];
wpt.pos.y   = [0 100 100 200 160];
wpt.time    = [0 40 60 80 100];

t = 0:1:max(wpt.time);                % time
x_p = pchip(wpt.time,wpt.pos.x,t);    % cubic Hermite inerpolation
y_p = pchip(wpt.time,wpt.pos.y,t);
x_s = spline(wpt.time,wpt.pos.x,t);   % spline interpolation
y_s = spline(wpt.time,wpt.pos.y,t);

% graphics
subplot(311)
plot(wpt.time,wpt.pos.x,'ro',t,x_p,'b','linewidth',2)
hold on, plot(t,x_s,'k'),hold off
grid, ylabel('x_d(t)')
legend('way-points','pchip','spline')

subplot(312)
plot(wpt.time,wpt.pos.y,'ro',t,y_p,'b','linewidth',2)
hold on, plot(t,y_s,'k'),hold off
grid, ylabel('y_d(t)')
legend('way-points','pchip','spline')

subplot(313)
plot(wpt.pos.y,wpt.pos.x,'ro',y_p,x_p,'b','linewidth',2)
hold on, plot(y_s,x_s,'k'),hold off
grid, xlabel('East (y_d)'),ylabel('North (x_d)')
legend('way-points','pchip','spline')
    
    