function SIMaidedINSheave()
% SIMaidedINSheave is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates an Inertial Navigation System (INS) aided by pressure 
% measurements using the Error-State Kalman Filter (ESKF). 
% 
% The position measurement frequency f_pos (typically 10 to 100 Hz) can be chosen 
% smaller or equal to the sampling frequency f_s (typically 1000 Hz), which is 
% equal to the Inertial Measurement Unit (IMU) measurement frequency. 
%
% Dependencies:
%   ins_heave.m     - Feedback ESKF for INS aided by pressure measurements. 
%   magneticField.m - Magnetic field vectors for different cities.
%  
% References:
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2024-11-02
% Revisions:

clear ins_heave; % Clear the persistent data structure 'ins'

%% USER INPUTS
T_final = 200;% Final simulation time (s)
f_s = 1000; % Sampling frequency equals IMU measurement frequency (Hz)
f_pos = 10; % Pressure measurement frequency (Hz)

% Sampling times
h  = 1/f_s; 	 
h_pos = 1/f_pos; 

% Constants
p_0 = 101325; % Air pressure in Pa at the surface
rho = 1025; % Desnity of water in kg/m^3
g = 9.81; % Acceleration of gravity in m/s^2

% Initialization of ESKF covariance matrix
P_prd = eye(3);
Rd = 0.001; % Pressure measurement covariance
Qd = diag([100 0.001]); % Velocity and acc bias covariances

% Initialization of INS states
z_ins = 0; 
v_z_ins = 0;
b_acc_ins = 0;
x_ins = [z_ins v_z_ins b_acc_ins];

% Initial values for signal generator
[m_ref, ~, mu,cityName] = magneticField(1); % Magntic field and latitude for city #1
displayMethod(cityName);
b_acc = [0.1 0.3 -0.1]'; % IMU biases
b_ars = [0.05 0.1 -0.05]';
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	% Initial states 

% Time vector initialization
t_slow = 0; % Initialize the time for the next slow measurement
t = 0:h:T_final; % Time vector from 0 to T_final          
nTimeSteps = length(t); % Number of time steps

%% MAIN LOOP
simdata = zeros(nTimeSteps,6); % Pre-allocate table for simulation data
posdata = zeros(floor(T_final * f_pos), 2); % Pre-allocate table for pos data
pos_index = 0; % Initialize index for posdata

for i=1:nTimeSteps
    
    % INS signal generator
    [x, f_imu, ~] = insSignal(x, h, t(i), mu, m_ref);
    phi = x(10); % roll angle
    theta = x(11); % pitch angle

    % Positions measurements are slower than the sampling time
    if t(i) > t_slow
        % Aiding
        pos_index = pos_index + 1; 
        posdata(pos_index, :) = [t(i), x(3)]; % Store position measurements 
        p = p_0 + rho * g * x(3) + 0.1 * randn; % Pressure measurement

        [x_ins, P_prd] = ins_heave(...
            x_ins, P_prd, h, Qd, Rd, f_imu, phi, theta, p_0, p);

        % Update the time for the next slow position measurement
        t_slow = t_slow + h_pos;  
    else  
        % No aiding
        [x_ins, P_prd] = ins_heave(x_ins, P_prd, h, Qd, Rd, f_imu, phi, theta);
    end

    % Store simulation data in a table (for testing)
    simdata(i,:) = [x(3) x(6) x(9) x_ins'];

end

%% PLOTS         
x     = simdata(:,1:3); % High-rate IMU data
x_hat = simdata(:,4:6); 

t_m = posdata(:,1); % Slow-rate position data
y_m = posdata(:,2);

figure(1); 

subplot(311)
h1 = plot(t_m,y_m,'xr'); hold on;
h2 = plot(t,x_hat(:,1),'b'); hold off;
xlabel('time (s)'),title('Down position [m]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_pos), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(312)
h1 = plot(t,x(:,2),'r'); hold on;
h2 = plot(t,x_hat(:,2),'b'); hold off;
xlabel('time (s)'),title('Down velocity [m/s]'),grid
legend([h1(1),h2(1)],['True down velocity at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(313)
h1 = plot(t,x(:,3),'r'); hold on;
h2 = plot(t,x_hat(:,3),'b'); hold off;
xlabel('time (s)'),title('Acc bias'),grid
legend([h1(1),h2(1)],['True acc bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)

%% DISPLAY DATA
function displayMethod(cityName)
    disp('-------------------------------------------------------------------');
    disp('MSS toolbox: Error-state Kalman filter (ESKF) for heave estimation');
    disp(['INS aided by pressure measurements at ',num2str(f_pos), ' Hz']);
    disp(['IMU measurements (specific force) at ',num2str(f_s),' Hz']);
    disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
    disp('-------------------------------------------------------------------');
    disp('Simulating...');
end

end