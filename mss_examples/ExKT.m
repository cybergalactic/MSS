% ExKT  Script for computation of Nomoto gain and time constants using 
%       nonlinear least-squares. The rudder input is 5 deg at t=0.
%
% Author:    Thor I. Fossen
% Date:      18th July 2001
% Revisions: 

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
