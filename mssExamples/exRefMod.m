% exRefMod   Second-order reference model with nonlinear damping and velocity 
%            saturation.
%
% Author:    Thor I. Fossen
% Date:      2001-11-03
% Revisions: 

zeta = 1;       % Relative damping ratio [-]
w = 1;          % Natural frequency [rad/s]
delta = 1;      % Quadratic damping coefficient [s/m]
v_max = 1;      % Max velocity [m/s]
h = 0.1;        % Sampling time [s]
T_final = 20;   % Simulation time [s]

% Time vector initialization
t = 0:h:T_final;               % Time vector from 0 to T_final          
nTimeSteps = length(t);        % Number of time steps

% Case 1: Linear mass-damper-spring system
r1 = 10 * ones(max(size(t)),1);
[A,B,C,D] = ord2(w,zeta); 
[x1,y1]   = lsim(A,B,C,D,r1,t); 

r2 = 10*r1;
[A,B,C,D] = ord2(w,zeta); 
[x2,y2]   = lsim(A,B,C,D,r2,t); 

% Case 2: Nonlinear damping
x = 0;
v = 0;
r2 = 10;
y2 = zeros(nTimeSteps,2);

for i=1:nTimeSteps
   y2(i,:) = [x v];   
   x_dot = v;
   v_dot = w^2 * (r2 - x) - 2 * zeta * w * v - delta *abs(v) * v;
   v = v + h * v_dot;
   x = x + h * x_dot;
end

% Case 3: Velocity saturation
x = 0;
v = 0;
r3 = 10;
y3 = zeros(nTimeSteps,2);

for i=1:nTimeSteps
   y3(i,:) = [x v];   

   v_dot = w^2 * (r3 -x) - 2 * zeta * w * v;
   x_dot = v;

   % Numerical integration
   v = v + h * v_dot;
   v = sat(v,v_max);  % Velocity saturation

   x = x + h * x_dot;
end

%% Plots
figure(gcf)
subplot(211); plot(t,y1(:,1))
hold on; plot(t,y2(:,1),'--k',t,y3(:,1),'-.r'); hold off; grid
title('Second-Order Reference Model with Nonlinear Damping and Velocity Saturation')
legend( ...
    'Linear damping', ...
    sprintf('Nonlinear damping (\\delta = %.1f)', delta), ...
    sprintf('Velocity saturation (v_{max} = %.1f m/s)', v_max), ...
    'Location','best')
subplot(212); plot(t,y1(:,2));
hold on; 
plot(t,y2(:,2),'--k',t,y3(:,2),'-.r'); 
yline(v_max,'-k','linewidth',1.0);
hold off; grid
legend( ...
    'Linear damping', ...
    sprintf('Nonlinear damping (\\delta = %.1f)', delta), ...
    sprintf('Velocity saturation (v_{max} = %.1f m/s)', v_max), ...
    'Location','best')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)


