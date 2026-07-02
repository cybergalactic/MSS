% exSpline is compatibel with MATLAB and GNU Octave (www.octave.org). 
% Cubic Hermite and spline interpolation of waypoints.
%
% Author:    Thor I. Fossen
% Date:      2002-07-06
% Revisions: 

% Waypoints 
wpt.pos.x   = [0 100 500 700 1000];
wpt.pos.y   = [0 100 100 200 160];
wpt.varpi   = [0 40 60 80 100];

varpi = 0:1:max(wpt.varpi);            % Path variable
x_p = pchip(wpt.varpi,wpt.pos.x,t);    % Cubic Hermite itnerpolation
y_p = pchip(wpt.varpi,wpt.pos.y,t);
x_s = spline(wpt.varpi,wpt.pos.x,t);   % Spline interpolation
y_s = spline(wpt.varpi,wpt.pos.y,t);

%% Plots
subplot(311)
plot(wpt.varpi,wpt.pos.x,'ro',varpi,x_p,'b')
hold on, plot(varpi,x_s,'k--'), hold off
grid, ylabel('x_d(\varpi)')
xlabel('Path parameter \varpi')
legend('Waypoints','pchip','spline')

subplot(312)
plot(wpt.varpi,wpt.pos.y,'ro',varpi,y_p,'b')
hold on, plot(varpi,y_s,'k--'), hold off
grid, ylabel('y_d(\varpi)')
xlabel('Path parameter \varpi')
legend('Waypoints','pchip','spline')

subplot(313)
plot(wpt.pos.y,wpt.pos.x,'ro',y_p,x_p,'b')
hold on, plot(y_s,x_s,'k--'),hold off
grid, xlabel('East (y_d)'),ylabel('North (x_d)')
legend('Waypoints','pchip','spline')
    
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)
set(findall(gcf,'type','line'),'linewidth',1.5)
    
