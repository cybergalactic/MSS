% exEKF is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates a discrete-time Extended Kalman filter (EKF), which
% demonstrates how the predictor-corrector representation can be applied to
% the nonlinear model:
%
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
%        y = x1 + white noise
%
% GNSS measurements are available at f_gnss [Hz], whereas the EKF predictor
% is propagated at the sampling frequency f_s [Hz]. The measurement update
% is applied only when a GNSS sample is available.
%
% Author:    Thor I. Fossen
% Date:      2018-10-17
% Revisions:

%% USER INPUTS
T_final = 100;  % Final simulation time (s)
f_s    = 10;    % Sampling frequency [Hz]
f_gnss = 1;     % GNSS position measurement frequency [Hz]

h      = 1/f_s;       % Sampling time (s)
h_gnss = 1/f_gnss;    % GNSS sampling time (s)

if f_gnss > f_s
    error("The GNSS frequency must be smaller than or equal to f_s.");
end

% Model parameters
%   dx1/dt = x2
%   dx2/dt = a * x2 * abs(x2) + b * u + white noise
a = -1;
b = 1;
e = 1;

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
ydata = [];                     % Measurement data: [time, measurement]
nextGnssTime = 0;               % Time of next GNSS measurement

for k = 1:nTimeSteps

   % Plant
   u = 0.1 * sin(0.1*t(k));
   w = 0.1 * randn;
   f = [ x(2)
         a * x(2) * abs(x(2)) + b * u ];
   E = [0 e]';
   x_dot = f + E * w;

   % Corrector: x_hat[k] and P_hat[k]
   newGNSSMeasurement = (t(k) >= nextGnssTime - 1e-10); 

   if newGNSSMeasurement

      y = x(1) + 0.1 * randn;
      ydata = [ydata; t(k), y];

      Cd = [1 0];

      S = Cd * P_prd * Cd' + Rd;
      K = P_prd * Cd' / S;
      IKC = eye(2) - K * Cd;

      x_hat = x_prd + K * (y - Cd * x_prd);
      P_hat = IKC * P_prd * IKC' + K * Rd * K';

      nextGnssTime = nextGnssTime + h_gnss;
   else
      x_hat = x_prd;     % No measurement update
      P_hat = P_prd;
   end

   % Store simulation data
   simdata(k,:) = [x' x_hat' P_hat(1,1) P_hat(2,2)];

   % Predictor: x_prd[k+1] and P_prd[k+1]
   f_hat = [ x_hat(2)
             a * x_hat(2) * abs(x_hat(2)) + b * u ];

   Ad = [ 1  h
          0  1 + h * 2 * a * abs(x_hat(2)) ];
   Ed = h * E;

   x_prd = x_hat + h * f_hat;
   P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

   % Euler's method: x[k+1]
   x = x + h * x_dot;
end

%% PLOTS
x     = simdata(:,1:2);
x_hat = simdata(:,3:4);
P_hat = simdata(:,5:6);

t_m = ydata(:,1);
y_m = ydata(:,2);

clf
figure(gcf)

subplot(211), plot(t_m,y_m,'xb',t,x_hat(:,1),'r')
xlabel('time (s)'), title('Position x_1'), grid
legend(['Measurement y = x_1 at ', num2str(f_gnss), ' Hz'], ...
       ['Estimate x_1hat at ', num2str(f_s), ' Hz']);

subplot(212), plot(t,x(:,2),'b',t,x_hat(:,2),'r')
xlabel('time (s)'), title('Velocity x_2'), grid
legend(['True velocity x_2 at ', num2str(f_s), ' Hz'], ...
       ['Estimate x_2hat at ', num2str(f_s), ' Hz']);

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)
set(findall(gcf,'type','line'),'linewidth',2)