% exEKFparameterEstimation is compatibel with MATLAB and GNU Octave
% (www.octave.org). This script simulates a discrete-time extended Kalman
% filter (EKF) implementation demonstrating how the "predictor-corrector
% representation" can be used to estimate unknown states and parameters
% in the nonlinear model:
% 
%   x1_dot = x2
%   x2_dot = a * x2 * abs(x2) + b * u + white noise
%   a_dot  = white noise
%        y = x1 + white noise
%
% where x3 = a is an unknown parameter.
%
% Author:    Thor I. Fossen
% Date:      14 July 2021
% Revisions: 

%% USER INPUTS
h  = 0.05;                      % sampling time [s]

% simulation parameters
N  = 4000;                      % no. of iterations
   
% model parameters 
a  = -0.9;                       % a is unknown in the EKF      
b  = 1;

% initial values for x and u
x = [ 0 0 ]';                    % x = [ x1 x2 ]
u = 0;

% initialization of the EKF
x_prd = [0 0 0]';                % x_prd = [ x1 x2 a ]'
P_prd = diag([1 1 1]);
Qd = diag([1 0.1]);               % Q = diag( Q_x2  Q_a )
Rd = 10 * diag (1);              % R = diag( R_x1 )

%% MAIN LOOP
simdata = zeros(N+1,7);                    % table of simulation data
                                           % simdata = [t x' x_hat' y' ]
for i=1:N+1
   t = (i-1) * h;                          % time (s)             

   % Nonlinear system
   u = 0.1 * sin(0.1*t);                   % input
   w = 0.1 * randn(1,1);                   % process noise
   f = [ x(2)
         a * x(2) * abs(x(2)) + b * u ];
   x_dot = f + [0 1]' * w;                 % dx/dt = f + E * w   
   
   % GNSS measurements 
   y = x(1) + 0.1 * randn(1);              % y = h + v
              
   % EKF measurement matrix Cd
   % where Cd = dh/dx is inearized about x = x_prd
   Cd     = [1 0 0]; 
   
   % KF gain       
   K = P_prd * Cd' * inv( Cd * P_prd * Cd' + Rd );
        
   % corrector   
   IKC = eye(3) - K * Cd;
   P_hat = IKC * P_prd * IKC' + K * Rd * K';
   x_hat = x_prd + K * ( y - Cd * x_prd );   
   
   % store simulation data in a table   
   simdata(i,:) = [t x' x_hat' y' ];    

   % discrete-time EKF model
   f_hat = [ x_hat(2)
             x_hat(3) * x_hat(2) * abs(x_hat(2)) + b * u 
             0 ];
      
   % Ad = I + h * A and Ed = h * E
   % where A = df/dx is linearized about x = x_hat
   Ad   = eye(3) + h * ...
       [ 0   1   0
         0   2 * x_hat(3) * abs(x_hat(2)) x_hat(2) * abs(x_hat(2)) 
         0   0   0 ];   
        
   % Ed = h * E; 
   Ed = h * [ 0 0
              1 0
              0 1 ];     
  
   % Predictor (k+1)       
   x_prd = x_hat + h * f_hat;
   P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';
   
   % Euler integration (k+1)
   x = x + h * x_dot;
end

%% PLOTS
t     = simdata(:,1); 
x     = simdata(:,2:3); 
x_hat = simdata(:,4:6); 
y_m   = simdata(:,7);

clf
figure(gcf)

subplot(311),plot(t,y_m,'b',t,x_hat(:,1),'r')
xlabel('time (s)'),title('Position x_1'),grid
legend('Measurement y = x_1','Estimate x_1hat');

subplot(312),plot(t,x(:,2),'b',t,x_hat(:,2),'r')
xlabel('time (s)'),title('Velocity x_2'),grid
legend('True veloicty x_2','Estimate x_2hat');

subplot(313),plot([t(1) t(end)],[a a],'b',t,x_hat(:,3),'r','linewidth',2)
xlabel('time (s)'),title('Parameter a'),grid
legend('True parameter a','Estimate x_3hat');

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
