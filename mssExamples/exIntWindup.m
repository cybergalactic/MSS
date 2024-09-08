%% exIntWindup
% exIntWindup is compatible with MATLAB and GNU Octave (www.octave.org).
% This script demonstrates integrator windup when the control law is
% saturated. The solution is to modify the integrator term using anti-windup. 
% The control objective is to regulate x to constant reference x_ref using:
%
% Continuous-time system:   
%                x_dot  = a * x + u + w
% Control law:            
%                e = x - x_ref
%                u = - Kp * e - Ki * z_dot,    z_dot = e
% 
% Definitions:   w: Disturbance 
%                u: Control input
%                x: State
%                z: Integral state
%
% Author: 2024-09-04 Thor I. Fossen

clearvars;

%% USER INPUTS
h  = 0.1;                   % Sample time in seconds
T_final = 200;              % Final simulation time in seconds

a  = -0.1;                  % Model parameter
            
Kp  = 1;                    % PI controller gains
Ki  = 0.5;
u_max = 0.7;                % Saturation 
u_min = -0.7;
t_switch = 100;             % Setpoint switching time

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

%% MAIN LOOP #1 (Without Anti-Windup)
table1 = zeros(nTimeSteps,3);  % Table for simulation data
x = 0; z = 0; x_ref = 5;

for i = 1:nTimeSteps
    
   % Disturbance
   w = 0.01 * randn(1);   

   % Switch setpoint
   if t(i) > t_switch
       x_ref = -2;
   end

   % PI control law 
   e = x - x_ref;
   u = -Kp * e - Ki * z;  
   
   % Saturation
   if u > u_max 
       u = u_max;  
   elseif u < u_min
       u = u_min;
   end
   
   % Store data in table1
   table1(i,:) = [x z u];                 
   
   % Euler integration   
   z = z + h * e;   
   x = x + h * (a * x + u + w);               
   
end
 
%% MAIN LOOP #2 (With Anti-Windup)
table2 = zeros(nTimeSteps,3);  % Table for simulation data
x = 0; z = 0; x_ref = 5; 

for i = 1:nTimeSteps
    
   % Disturbance
   w = 0.01 * randn(1);  

   % Switch setpoint
   if t(i) > t_switch
       x_ref = -2;
   end

   % Unsaturated PI control law 
   e = x - x_ref;
   u_unsat = -Kp * e - Ki * z;              
   
   % No anti-wind up
   u = u_unsat;
   z = z + h * e;

   % Saturated control law and integrator anti-windup
   if u_unsat > u_max
       u = u_max; % Saturation
       z = z - (h/Ki) * (u - u_unsat); % Anti-windup
   elseif u < u_min
       u = u_min; % Saturation
       z = z - (h/Ki) * (u - u_unsat); % Anti-windup
   end
   
   % Store data in table2
   table2(i,:) = [x z u];                 
   
   % Euler integration    
   x = x + h * (a * x + u + w);               
   
end

%% PLOTS  
x1   = table1(:,1);  
z1   = table1(:,2); 
u1   = table1(:,3); 
 
x2   = table2(:,1);  
z2   = table2(:,2); 
u2   = table2(:,3);  

figure(gcf)
subplot(311),plot(t,x1,t,x2,'linewidth',2),xlabel('time (s)'),title('state x'),grid
legend('without anti-wind-up','with anti-wind-up')
subplot(312),plot(t,z1,t,z2,'linewidth',2),xlabel('time (s)'),title('integral state z'),grid
legend('Without anti-windup','With anti-windup')
subplot(313),plot(t,u1,t,u2,'linewidth',2),xlabel('time (s)'),title('control u'),grid
legend('Without anti-windup','With anti-windup')

