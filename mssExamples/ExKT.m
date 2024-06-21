% exKT requires the MATLAB Optimization Toolbox.
% Script for computation of the Nomoto gain and time constants using 
% nonlinear least-squares. The rudder input is 5 deg at t = 0.
%
% Author:    Thor I. Fossen
% Date:      2001-07-18
% Revisions: 

N = 2000;  % Number of samples
h = 0.1;   % Sample time

xout = zeros(N,2); 
x = zeros(7,1);            % x = [ u v r psi x y delta ]' 
delta_R = 5*(pi/180);      % Rudder angle step input

for i=1:N
    xout(i,:) = [(i-1)*h ,x(3)]; 
    xdot = mariner(x,delta_R);   % Nonlinear Mariner model
    x = euler2(xdot,x,h);        % Euler integration
end

% Time-series
tdata = xout(:,1);
rdata = xout(:,2)*180/pi;   

% Nonlinear least-squares parametrization: T dr/dt + r = K delta, delta = -delta_R
% x(1) = 1/T and x(2) = K
x0 = [0.1 1]'
F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*5','x','tdata')
x = lsqcurvefit(F,x0, tdata, rdata);

% Estimated parameters
T = 1/x(1)
K = x(2)

figure(gcf),subplot(211)
plot(tdata,rdata,'g',tdata, ...
    exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*5, ...
    'r', 'linewidth',2)
grid
title('Nonlinear least-squares fit of Mariner model for \delta = 5 (deg)'),xlabel('time (s)')
legend('Nonlinear model','Estimated 1st-order Nomoto model')
