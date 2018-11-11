% ExKF Discrete-time Kalman filter (KF) implementation demonstrating
% how the "predictor-corrector representation" can be applied to the
% following linear model:
% 
%   dx/dt = A * x + B * u + E * w 
%       y = H * x + v
%
%   x_k+1 = Ad * x_k + Bd * u_k + Ed * w_k 
%     y_k = H * x_k + v_k
%
% Let h be the sampling time and Ad = I + h * A, Bd = h * B and Ed = h * E. 
%
% The case study is a mass-damper system with position measured at frequency
% f_pos [Hz], which can be chosen smaller or equal to the sampling frequency 
% f_samp [Hz]. The ratio between the frequencies must be an integer:
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

% model paramters for mass-damper system (x1 = position, x2 = velcoity)
A = [0 1
     0 -0.1 ];
B = [0
     1];
E = [0
     1];
H = [1 0]; 

% discrete-time matrices
Ad = eye(2) + h * A;
Bd = h * B;
Ed = h * E;
 
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
   x_dot = A * x + B * u + E * w;          
     
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
      
   % Predictor (k+1)  
   x_prd = Ad * x_hat + Bd * u;
   P_prd = Ad * P_hat * Ad' + Ed * Q * Ed';
   
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