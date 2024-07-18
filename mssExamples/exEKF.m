% exEKF is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates a discrete-time Extended Kalman filter (EKF), which 
% demonstrates how the "predictor-corrector representation" can be 
% applied to the nonlinear model:
% 
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
%        y = x1 + white noise
%
% The GNSS position measurement frequency f_gnss [Hz] can be chosen smaller 
% or equal to the  sampling frequency f_s [Hz]. The ratio between the 
% frequencies must be an integer:
%
%     Integer:  Z = f_s/f_gnss >= 1 
%
% Author:    Thor I. Fossen
% Date:      2018-10-17
% Revisions: 

%% USER INPUTS
T_final = 100;  % Final simulation time (s)
f_s    = 10;    % sampling frequency [Hz]
f_gnss = 1;     % GNSS position measurement frequency [Hz]

Z = f_s/f_gnss;
if ( mod(Z,1) ~= 0 || Z < 1 )
    error("f_s = Z * f_m is not specified such that Z = f_s/f_m >= 1%"); 
end

h  = 1/f_s; 	  % sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss;      

% Model parameters 
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
a  = -1;              
b  = 1;
e  = 1;

% Initial values for x and u
x = [0 0]';	        
u = 0;

% Initialization of the Kalman filter
x_prd = [0 0]';        
P_prd = diag([1 1]);
Qd = 1;
Rd = 10;

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

%% MAIN LOOP
simdata = zeros(nTimeSteps,6);  % Pre-allocate table of simulation data
ydata = [0 x_prd(1)];           % Pre-allocate table of measurement data

for i=1:nTimeSteps
           
   % Plant
   u = 0.1 * sin(0.1*t(i));                % Input
   w = 0.1 * randn(1);                     % Process noise
   f = [ x(2)
         a*x(2)*abs(x(2)) + b*u];
   E = [0 e]';
   x_dot = f + E * w;                      % dx/dt = f + E * w   
   
   % GNSS measurements are Z times slower than the sampling time
   if mod( t(i), h_gnss ) == 0
       y = x(1) + 0.1 * randn(1); 
       ydata = [ydata; t(i), y]; 
       
       Cd     = [1 0];               
   else
       Cd     = [0 0];               % No measurement
   end
    
   % KF gain      
   K = P_prd * Cd' * inv( Cd * P_prd * Cd' + Rd );
        
   % Corrector   
   IKC = eye(2) - K*Cd;
   P_hat = IKC * P_prd * IKC' + K * Rd * K';
   eps = y - Cd * x_prd;
   x_hat = x_prd + K * eps;   
   
   % Store simulation data in a table   
   simdata(i,:) = [x' x_hat' P_hat(1,1) P_hat(2,2) ];    

   % Discrete-time EKF model
   f_hat = [ x_hat(2)
             a * x_hat(2) * abs(x_hat(2)) + b * u ];
   f_d   = x_hat + h * f_hat;
      
   % Predictor (k+1)  
   % Ad = I + h * A and Ed = h * E
   % where A = df/dx is linearized about x = x_hat
   Ad   = [ 1   h
             0  1 + h*2*a*abs(x_hat(2)) ];   
   Ed = h * E; 
   
   x_prd = f_d;
   P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';
   
   % Euler's method (k+1)
   x = x + h * x_dot;
end

%% PLOTS
x     = simdata(:,1:2); 
x_hat = simdata(:,3:4); 
X_hat = simdata(:,5:6);

t_m = ydata(:,1);
y_m = ydata(:,2);

clf
figure(gcf)

subplot(211),plot(t_m,y_m,'xb',t,x_hat(:,1),'r')
xlabel('time (s)'),title('Position x_1'),grid
legend(['Measurement y = x_1 at ', num2str(f_gnss), ' Hz'],...
    ['Estimate x_1hat at ', num2str(f_s), ' Hz']);

subplot(212),plot(t,x(:,2),'b',t,x_hat(:,2),'r')
xlabel('time (s)'),title('Velocity x_2'),grid
legend(['True veloicty x_2 at ', num2str(f_s), ' Hz'],...
    ['Estimate x_2hat at ', num2str(f_s), ' Hz']);

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
