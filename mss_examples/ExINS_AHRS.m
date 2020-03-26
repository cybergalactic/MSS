% ExINS_AHRS error-state (indirect) feedback Kalman filter for INS aided by
% position measurements. It is assumed that the attitude (phi, theta, psi) is 
% measured by an AHRS. Alternatively, the roll and pitch angles can be 
% computed from specific force measurements (to less accuracy) using: 
%
%    [phi, theta] = acc2rollpitch(f)
%
% while the yaw angle, psi, is measured by a compass.
%
% The GNSS position measurement frequency f_gnss can be chosen smaller or
% equal to the  sampling frequency f_s, which is equal to the IMU
% meaurement frequency. The ratio between the frequencies must be an integer:
%
%     Integer:  Z = f_s/f_gnss >= 1 
%
% The main loop calls:
%
% [x_ins, P_prd] = ins_ahrs( ...
%      x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_ahrs, y_pos, y_vel)
%
% each time a GNSS position y_pos is received at the slow frequency f_gnss.
% For samples without new measurements, the arguments y_pos and y_vel are
% omitted when calling the function ins_ahrs when there are no new GNSS 
% measurements.
%
% Author:    Thor I. Fossen
% Date:      21 March 2020
% Revisions: 26 March 2020, added flags as user input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_s    = 100;   % sampling frequency [Hz]
f_gnss = 1;     % GNSS measurement frequency [Hz]

% Flag
vel = 0;          % 0 = no velocity meaurement, 1 = velocity aiding
compass = 0;      % 0 = AHRS, 1 = compass and [phi, theta] = acc2roll(f_imu)

% Parameters
Z = f_s/f_gnss;   % ratio betwween sampling/IMU frequencies
h  = 1/f_s; 	  % sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss;  

% simulation parameters
N  = 10000;		  % no. of iterations
b_acc = [0.1 0.3 -0.1]';
b_ars = [0.05 0.1 -0.05]';
   
% initial values for x (for testing)
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	        

% initialization of Kalman filter
P_prd = eye(15);

if (vel == 0)
   Rd = diag([1 1 1  1 1 1]);       % only position measurements
   Qd = diag([1 1 1  1 1 1  10 10 10  0.01 0.01 0.01]);
else
   Rd = diag([10 10 10 1 1 1 0.1 0.1 0.1]);  % velocity measurements
   Qd = diag([1 1 1  1 1 1  0.1 0.1 0.1  0.01 0.01 0.01]);
end

% initialization of INS
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
th_ins = [0 0 0]';
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; th_ins; b_ars_ins];

% WGS-84 gravity model
mu = 63.4305 * pi / 180;    % lattitude  
g = gravity(mu);  


%% Display
disp('----------------------------------------------------------');
disp('MSS toolbox: Error-state (indirect) feedback Kalman filter');
disp('Attitude parametrization: Euler angles');
if (vel == 0)
   disp(['INS aided by position at ',num2str(f_gnss), ' Hz']);
else
    disp(['INS aided by position and velocity at ',num2str(f_gnss), ' Hz']);  
end
   disp(['IMU measurements (specific force and ARS) at ',num2str(f_s),' Hz']);
if (compass == 0)
   disp(['AHRS measurements (phi, theta, psi) at ',num2str(f_s), ' Hz']);
else
   disp(['COMPASS measurements (psi) at ',num2str(f_s), ' Hz']);
end
disp('----------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(N+1,34);                  % table of simulation data
ydata = [0 x(1:3)'];                      % table of position measurements

for i=1:N+1
    t = (i-1) * h;                        % time (s)   
    
    % Signal generator (for testing)
    g_n = [0 0 g]';                       % NED gravity vector
    g_b = Rzyx(x(10),x(11),x(12))' * g_n; % BODY gravity vector
    
    f_true = [0.1  * sin(0.1*t)           % true specific force: f = a - g 
              0.1  * cos(0.1*t)
              0.05 * sin(0.05*t)] - g_b;

    w_true = [ 0.01 * cos(0.2*t)           % true angular rate
              -0.05 * sin(0.1*t)
               0.05 * sin(0.1*t) ];
       
    x_dot = [ x(4:6)                                   % true states
              Rzyx(x(10),x(11),x(12)) * f_true + g_n 
              zeros(3,1)    
              Tzyx(x(10),x(11)) * w_true 
              zeros(3,1)  ];
   
    % IMU measurements          
    w1 = 0.01 * randn(3,1);                 % white noise     
    w2 = 0.01 * randn(3,1);    
    f_imu = f_true + b_acc + w1;   
    w_imu = w_true + b_ars + w2;
    
    % AHRS measurements
    y_theta = x(10:12) + 0.01 * randn(3,1); % true roll, pitch, yaw  
    y_ahrs = y_theta;
    
    % Use specific force and compass instead of an AHRS
    if (compass == 1)
       [phi, theta] = acc2rollpitch( f_imu );
        y_ahrs = [phi, theta, x(12)]';
    end
    
    % GNSS measurements are Z times slower than the sampling time
    if mod( t, h_gnss ) == 0
        
        y_pos = x(1:3) + 0.1 * randn(3,1);      % position measurements
        y_vel = x(4:6) + 0.01 * randn(3,1);     % optionally velocity meas.
        ydata = [ydata; t, y_pos'];             % store position measurements                  
        
        if (vel == 0)
            [x_ins,P_prd] = ins_ahrs(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs,y_pos);
        else
            [x_ins,P_prd] = ins_ahrs(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs,y_pos,y_vel);
        end
        
    else  % no aiding
        
        [x_ins,P_prd] = ins_ahrs(x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs);
        
    end
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x' x_ins', y_theta']; 
    
    % Signal generator (for testing)
    x = x + h * x_dot;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t     = simdata(:,1);           
x     = simdata(:,2:16); 
x_hat = simdata(:,17:31); 

theta_m = simdata(:,32:34);     % fast AHRS measurements

t_m = ydata(:,1);               % slow GNSS measurements
y_m = ydata(:,2:4);

figure(1); figure(gcf)

subplot(311),plot(t_m,y_m,'xb',t,x_hat(:,1:3),'r')
xlabel('time (s)'),title('Position [m]'),grid
legend(['Measurement at ', num2str(f_gnss), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(312),plot(t,x(:,4:6),'b',t,x_hat(:,4:6),'r')
xlabel('time (s)'),title('Velocity [m/s]'),grid
legend(['True velocity at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(313),plot(t,x(:,7:9),'b',t,x_hat(:,7:9),'r')
xlabel('time (s)'),title('Acc bias'),grid
legend(['True acc bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

figure(2); figure(gcf)

subplot(211),plot(t,(180/pi)*theta_m,'b',t,(180/pi)*x_hat(:,10:12),'r')
xlabel('time (s)'),title('Angle [deg]'),grid
legend(['Measurement at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(212),plot(t,x(:,13:15),'b',t,x_hat(:,13:15),'r')
xlabel('time (s)'),title('ARS bias'),grid
legend(['True ARS bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

