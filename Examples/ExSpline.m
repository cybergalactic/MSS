% ExSpline   Cubic Hermite and spline interpolation of way-points 
% Author:    Thor I. Fossen
% Date:      6 July 2002
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org

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
    
    