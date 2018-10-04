% ExKT  Script for computation of Nomoto gain and time constants using 
%       nonlinear least-squares. The rudder input is 5 deg at t=0.
%
% Author:    Thor I. Fossen
% Date:      18th July 2001
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
%
N = 2000;  % number of samples
h = 0.1;   % sample time

xout = zeros(N,2); 
x = zeros(7,1);            % x = [ u v r psi x y delta ]' 
delta_R = 5*(pi/180);      % rudder angle step input

for i=1:N,
    xout(i,:) = [(i-1)*h ,x(3)]; 
    xdot = mariner(x,delta_R);   % nonlinear Mariner model
    x = euler2(xdot,x,h);        % Euler integration
end

% time-series
tdata = xout(:,1);
rdata = xout(:,2)*180/pi;   

% nonlinear least-squares parametrization: T dr/dt + r = K delta,   delta = -delta_R
% x(1) = 1/T and x(2) = K
x0 = [0.1 1]'
F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*5','x','tdata')
x = lsqcurvefit(F,x0, tdata, rdata);

% estimated parameters
T = 1/x(1)
K = x(2)

figure(gcf),subplot(211)
plot(tdata,rdata,'g',tdata,exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*5,'r'),grid
title('Nonlinear least-squares fit of Mariner model for \delta = 5 (deg)'),xlabel('time (s)')
legend('Nonlinear model','Estimated 1st-order Nomoto model')
