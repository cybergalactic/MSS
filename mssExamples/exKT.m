% exKT requires the MATLAB Optimization Toolbox.
% Script for computation of the Nomoto gain and time constants using 
% nonlinear least-squares. The differential equation 
%
%   r_dot + (1/T) * r = (K/T) * delta 
%
% where delta = -delta_R is a step input, has the solution
% 
%   r(t) = K * delta_R * ( 1 - e^(-(t/T) )  
%
% Author:    Thor I. Fossen
% Date:      2001-07-18
% Revisions: 
%   2024-06-21 - Updated to avoid inline functions

N = 2000;  % Number of samples
h = 0.1;   % Sample time

xout = zeros(N,2); 
x = zeros(7,1);     % x = [ u v r psi x y delta ]' 
delta_R = 5;        % Rudder angle step input

for i = 1:N
    xout(i,:) = [(i-1)*h ,x(3)]; 
    xdot = mariner(x,deg2rad(delta_R));   % Nonlinear Mariner-class model
    x = euler2(xdot, x, h);               % Euler integration
end

% Time-series
tdata = xout(:,1);
rdata = rad2deg(xout(:,2)); 

% Nonlinear least-squares: T r_dot + r = K delta, delta = -delta_R
% x(1) = 1/T and x(2) = K
x0 = [0.1 1];

F = @(x, tdata) x(2) * (1 - exp(-tdata * x(1))) * delta_R;

x = lsqcurvefit(F, x0, tdata, rdata);

% Estimated parameters
T = 1 / x(1);
K = x(2);

% Display the results
disp(['Estimated K: ', num2str(K)]);
disp(['Estimated T: ', num2str(T)]);

% Plot the results
figure;
subplot(2,1,1);
plot(tdata, rdata, 'k', 'LineWidth', 2); hold on;
plot(tdata, x(2) * (1 - exp(-tdata * x(1))) * delta_R, 'r', 'LineWidth', 2);
grid on;
title('Nonlinear least-squares fit of Mariner-class vessel model for \delta = 5 deg');
xlabel('Time (s)');
ylabel('r (deg/s)');
legend('Nonlinear model', 'Estimated 1st-order Nomoto model','Location','southeast');

subplot(2,1,2);
plot(tdata, rdata - (x(2) * (1 - exp(-tdata * x(1))) * delta_R), 'LineWidth', 2);
grid on;
title('Residuals');
xlabel('Time (s)');
ylabel('Residuals (deg/s)');

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)