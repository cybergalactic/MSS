% exINS_AHRS is compatible with MATLAB and GNU Octave (www.octave.org).
% Error-state (indirect) feedback Kalman filter for INS aided by
% position measurements. It is assumed that the attitude (phi, theta, psi)
% is measured by an Attitude Heading Reference System (AHRS). Alternatively,
% the roll and pitch angles can be computed from specific force 
% measurements (to less accuracy) using: 
%
%    [phi, theta] = acc2rollpitch(f)
%
% while the yaw angle, psi, is measured by a compass. A more accurate and
% recommended solution is to use the MEKF; see exINS_MEKF. 
%
% The GNSS position measurement frequency f_gnss can be chosen smaller or
% equal to the  sampling frequency f_s, which is equal to the IMU
% measurement frequency. The ratio between the frequencies must be an integer:
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
%            28 March 2020, modified to use an INS signal generator

%% USER INPUTS
T_final = 100;	% Final simulation time (s)
f_s    = 100;   % Sampling frequency [Hz]
f_gnss = 1;     % GNSS measurement frequency [Hz]

% Flag
vel = 0;          % 0 = no velocity meaurement, 1 = velocity aiding
compass = 0;      % 0 = AHRS, 1 = compass and [phi, theta] = acc2roll(f_imu)

% Parameters
Z = f_s/f_gnss;   % ratio betwween sampling/IMU frequencies
h  = 1/f_s; 	  % sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss;  

% Biases
b_acc = [0.1 0.3 -0.1]';
b_ars = [0.05 0.1 -0.05]';
   
% Initial values for x (for testing)
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	        

% Initialization of Kalman filter
P_prd = eye(15);

if (vel == 0)
   Rd = diag([1 1 1  1 1 1]);       % only position measurements
   Qd = diag([1 1 1  1 1 1  10 10 10  0.01 0.01 0.01]);
else
   Rd = diag([10 10 10 1 1 1 0.1 0.1 0.1]);  % velocity measurements
   Qd = diag([1 1 1  1 1 1  0.1 0.1 0.1  0.01 0.01 0.01]);
end

% Initialization of INS
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
th_ins = [0 0 0]';
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; th_ins; b_ars_ins];

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

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

%% MAIN LOOP
simdata = zeros(nTimeSteps,33); % Pre-allocate table of simulation data
ydata = [0 x(1:3)'];            % Pre-allocatetable of position measurements

for i=1:nTimeSteps
    
    % INS signal generator  
    [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t(i));
    
    if (compass == 0)      % AHRS measurements
       
       y_ahrs = x(10:12);  % True roll, pitch, yaw                  %
       
    else                % Use specific force and compass instead of an AHRS
        
       [phi, theta] = acc2rollpitch( f_imu );
        y_ahrs = [ phi, theta, x(12) ]';
        
    end
    
    % GNSS measurements are Z times slower than the sampling time
    if mod( t(i), h_gnss ) == 0
        
        y_pos = x(1:3) + 0.1 * randn(3,1);    % position measurements
        y_vel = x(4:6) + 0.01 * randn(3,1);   % optionally velocity meas.
        ydata = [ydata; t(i), y_pos'];           % store position measurements                  
       
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
    simdata(i,:) = [x' x_ins', y_ahrs']; 
    
end

%% PLOTS       
x     = simdata(:,1:15); 
x_hat = simdata(:,16:30); 

theta_m = simdata(:,31:33);     % fast AHRS measurements

t_m = ydata(:,1);               % slow GNSS measurements
y_m = ydata(:,2:4);

figure(1); figure(gcf)

subplot(311)
h1 = plot(t_m,y_m,'xr'); hold on;
h2 = plot(t,x_hat(:,1:3),'b'); hold off;
xlabel('time (s)'),title('Position [m]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_gnss), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz']);

subplot(312)
h1 = plot(t,x(:,4:6),'r'); hold on;
h2 = plot(t,x_hat(:,4:6),'b'); hold off;
xlabel('time (s)'),title('Velocity [m/s]'),grid
legend([h1(1),h2(1)],['True velocity at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(313)
h1 = plot(t,x(:,7:9),'r'); hold on;
h2 = plot(t,x_hat(:,7:9),'b'); hold off;
xlabel('time (s)'),title('Acc bias'),grid
legend([h1(1),h2(1)],['True acc bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2); figure(gcf)

subplot(211)
h1 = plot(t,rad2deg(theta_m),'b'); hold on;
h2 = plot(t,rad2deg(x_hat(:,10:12)),'r'); hold off;
xlabel('time (s)'),title('Angle [deg]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(212)
h1 = plot(t,x(:,13:15),'b'); hold on;
h2 = plot(t,x_hat(:,13:15),'r'); hold off;
xlabel('time (s)'),title('ARS bias'),grid
legend([h1(1),h2(1)],['True ARS bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

