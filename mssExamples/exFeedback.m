% exFeedback is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates numerical integration of a 1st-order system with 
% feedback and feedforward control laws. Euler's method is implemented 
% using a for-end loop, and the results are stored in a table and plotted.
%
% Continious-time system:    .
%                            x = a * x + b * u + w
% Control law:
%                            u = -k * x + kr * r
% 
% Definitions:               w = white noise
%                            u = control input
%                            x = state
%                            r = reference signal
%
% Author:                    2018-10-21 Thor I. Fossen

%% USER INPUTS
h  = 0.1;               % sample time (s)
N  = 500;               % number of samples

a  = -0.1;              % model parameters
b  = 1.0;

k  = 0.1;               % feedback gain
kr = -(a-b*k)/b;        % steady-state feedforward gain (x_ss = r)

r  = 0.5;               % reference signal (user input)

x = -0.1;               % initial state x(0)

table = zeros(N+1,3);   % table for simulation data

%% FOR LOOP
for i = 1:N+1
   t = (i-1)*h;             % time
   u = -k*x + kr*r;         % control law
   w = 0.02 * randn(1);     % process noise
   x_dot = a*x + b*u + w;   % system model
   
   table(i,:) = [t x u];    % store data in table
   
   x = x + h * x_dot;	    % Euler integration
end

%% PLOT FIGURES
t   = table(:,1);  
x   = table(:,2);  
u   = table(:,3);  

figure(gcf)
subplot(211),plot(t,x,'linewidth',2),xlabel('time (s)'),title('x'),grid
hold on
plot(t,r*ones(N+1,1),'linewidth',2)
hold off
subplot(212),plot(t,u,'linewidth',2),xlabel('time (s)'),title('u'),grid

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
