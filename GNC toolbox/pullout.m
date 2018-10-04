function [t,r1,r2] = pullout(ship,x,ui,h)
% PULLOUT  [t,r1,r2] = pullout(ship,x,ui,h)
%          pullout maneuver for a ship, see ExPullout.m
%
% Inputs :
% 'ship'          = ship model. Compatible with the models under .../gnc/VesselModels/
% x               = initial state vector for ship model
% ui              = [delta,:] where delta is the rudder command at time = t_rudderexecute
% h               = sampling time
%
% Outputs :
% t               = time vector
% r1,r2 =         = yaw rates for positiv and negative pullouts
%
% Author:   Thor I. Fossen
% Date:     25th July 2001
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

if nargin~=4, error('number of inputs must be 4'); end

T = 600;                 % maximum time for steady-state yaw rate to be reached
N = round((2*T)/h);      % number of samples corresponding to t=2T (sec)
xout1 = zeros(N+1,8);    % memory allocation
xout2 = zeros(N+1,8);    % memory allocation
xinit = x;
delta_c = ui(1);

disp('Simulating...')

% positive rudder step
u_ship = ui;

for i=1:N+1,
    time = (i-1)*h;
   
    if round(time) < T, 
        u_ship(1) = delta_c;
    else
        u_ship(1) = 0;
    end     
    
    [xdot,U] = feval(ship,x,u_ship);       % ship model
    
    xout1(i,:) = [time,x(1:6)',U];  
    
    x = euler2(xdot,x,h);                     % Euler integration
end

% negative rudder step
x = xinit;

for i=1:N+1,
    time = (i-1)*h;
   
    if round(time) < T, 
        u_ship(1) = -delta_c;
    else
        u_ship(1) = 0;
    end     
    
    [xdot,U] = feval(ship,x,u_ship);       % ship model
    
    xout2(i,:) = [time,x(1:6)',U];  
    
    x = euler2(xdot,x,h);                  % Euler integration
end

% time-series
t     = xout1(:,1);
r1    = xout1(:,4)*180/pi; 
r2    = xout2(:,4)*180/pi; 

% plots
figure(1)
subplot(111),
plot(t,r1,'b','linewidth',2)
hold on
plot(t,r2,'b','linewidth',2)
plot(t,0*zeros(length(t),1),'r','linewidth',2)
hold off
xlabel('time (s)'),title('yaw rate r (deg/s)'),grid


