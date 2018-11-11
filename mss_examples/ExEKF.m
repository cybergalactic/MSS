% ExEKF Discrete-time extended Kalman filter (EKF) implementation demonstrating
% how the "predictor-corrector representation" can be applied to the
% following nonlinear model:
% 
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
%        y = x1 + white noise
%
% The position measurement frequency f_pos [Hz] can be chosen smaller or
% equal to the  sampling frequency f_samp [Hz]. The ratio between the 
% frequencies must be an integer:
%
%      N_integer = f_samp/f_pos >= 1 
%
% Author:    Thor I. Fossen
% Date:      17 October 2018
% Revisions: 

%% USER INPUTS
f_samp = 10;    % sampling frequency [Hz]
f_pos = 0.5;    % GNSS position measurement frequency [Hz]

N_integer = f_samp/f_pos;
if ( mod(N_integer,1) ~= 0 || N_integer < 1 )
    error("f_samp = N_integer * f_pos is not specified such that N_integer >= 1"); 
end

% simulation parameters
N  = 1000;		  % no. of iterations

h  = 1/f_samp; 	  % sampling time: h  = 1/f_samp (s) 
h_pos = 1/f_pos;

% model parameters 
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
a  = -1;              
b  = 1;
e  = 1;

% initial values for x and u
x = [0 0]';	        
u = 0;

% initialization of Kalman filter
x_prd = [0 0]';        
P_prd = diag([1 1]);
Q = 1;
R = 10;

N_integer = f_samp/f_pos;
if ( mod(N_integer,1) ~= 0 || N_integer < 1 )
    error("f_samp = N_integer * f_pos is not specified such that N_integer >= 1"); 
end

%% MAIN LOOP
simdata = zeros(N+1,7);                    % table for simulation data
posdata = [0 x_prd(1)];                    % table for position measurements

for i=1:N+1
   t = (i-1) * h;                          % time (s)             

   % Plant
   u = 0.1 * sin(0.1*t);                   % input
   w = 0.1 * randn(1);                     % process noise
   f = [ x(2)
         a*x(2)*abs(x(2)) + b*u];
   E = [0 e]';
   x_dot = f + E * w;                      % dx/dt = f + E * w
   
   % Discrete-time linearized matrices: PHI = I + h * A and GAMMA = h * B 
   % where A = df/dx is linearized about x = x_prd
   PHI   = [ 1   h
             0  1 + h*2*a*abs(x_prd(2)) ];   
   GAMMA = h * E;    
   
   % measurements are slower than the sampling time
   if mod( t, h_pos ) == 0
       z = x(1) + 0.1 * randn(1);  
       posdata = [posdata; t z];    % store measurement in a table
       H     = [1 0];               
   else
       H     = [0 0];               % no measurement
   end
    
   % KF gain      
   K = P_prd * H' * inv( H * P_prd * H' + R );
        
   % corrector   
   IKH = eye(2) - K*H;
   P_hat = IKH * P_prd * IKH' + K * R * K';
   eps = z - H * x_prd;
   x_hat = x_prd + K * eps;   
   
   % store simulation data in a table   
   simdata(i,:) = [t x' x_hat' P_hat(1,1) P_hat(2,2) ];    

   % discrete-time extended KF-model
   f_hat = [ x_hat(2)
             a * x_hat(2) * abs(x_hat(2)) + b * u ];
   f_k   = x_hat + h * f_hat;
      
   % Predictor (k+1)  
   x_prd = f_k;
   P_prd = PHI * P_hat * PHI' + GAMMA * Q * GAMMA';
   
   % Euler integration (k+1)
   x = x + h * x_dot;
end

%% PLOTS
t     = simdata(:,1); 
x     = simdata(:,2:3); 
x_hat = simdata(:,4:5); 
X_hat = simdata(:,6:7);

t_pos = posdata(:,1);
z_pos = posdata(:,2);

clf
figure(gcf)

subplot(211),plot(t_pos,z_pos,'xb',t,x_hat(:,1),'r')
xlabel('time (s)'),title('Position x_1'),grid
legend(['Measurement z = x_1 at ', num2str(f_pos), ' Hz'],...
    ['Estimate x_1hat at ', num2str(f_samp), ' Hz']);

subplot(212),plot(t,x(:,2),'b',t,x_hat(:,2),'r')
xlabel('time (s)'),title('Velocity x_2'),grid
legend(['True veloicty x_2 at ', num2str(f_samp), ' Hz'],...
    ['Estimate x_2hat at ', num2str(f_samp), ' Hz']);