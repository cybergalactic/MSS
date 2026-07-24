function SIMaidedINSheave()
% SIMaidedINSheave is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates an Inertial Navigation System (INS) aided by pressure 
% measurements: 
% 
%   p = p_0 + rho * g * z
%
% using the Error-State Kalman Filter (ESKF). The pressure measurement frequency 
% f_pos (typically 10 to 100 Hz) can be chosen smaller or equal to the sampling 
% frequency f_s (typically 1000 Hz), which is equal to the Inertial Measurement 
% Unit (IMU) measurement frequency. 
%
% Dependencies:
%   ins_heave.m     - Feedback ESKF for INS aided by pressure measurements. 
%   magneticField.m - Magnetic field vectors for different cities.
%  
% References:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 3rd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2024-11-02
% Revisions:

clear ins_heave; % Clear persistent data structure 'ins' in the function

% ==============================================================================
% Simulation parameters
% ==============================================================================
T_final = 200; % Final simulation time (s)
f_s = 1000; % Sampling frequency equals IMU measurement frequency (Hz)
f_pos = 10; % Pressure measurement frequency (Hz)

% Sampling times
h  = 1/f_s; 	 
h_pos = 1/f_pos; 

% Constants
p_0 = 101325; % Air pressure in Pa at the surface
rho = 1025; % Density of water in kg/m^3
[m_ref, ~, mu,cityName] = magneticField(1); % Magnetic field and latitude for city #1
g = gravity(mu); % Acceleration of gravity in m/s^2

displayMethod(cityName, f_s, f_pos);

% Initialization of ESKF
T_acc = 300;                         % Accelerometer-bias time constant [s]

% Pressure/depth measurement noise
sigma_z = 0.01;                      % Depth standard deviation [m]
sigma_p = rho * g * sigma_z;         % Pressure standard deviation [Pa]
Rd = sigma_z^2;                      % Depth measurement variance [m^2]

% Initial error covariance
sigma_z0 = 0.5;                      % Initial depth uncertainty [m]
sigma_v0 = 0.1;                      % Initial vertical-velocity uncertainty [m/s]
sigma_b0 = 0.05;                     % Initial accelerometer-bias uncertainty [m/s^2]

P_prd = diag([sigma_z0^2, sigma_v0^2, sigma_b0^2]);

% Process-noise standard deviations
sigma_f = 0.1;                       % Specific-force noise [m/s^2]
sigma_b = 1e-4;                      % Bias-driving noise [m/s^2]

Qd = diag([sigma_f^2, sigma_b^2]);

% Initialization of INS states
z_ins = 0; 
v_z_ins = 0;
b_acc_ins = 0;
x_ins = [z_ins v_z_ins b_acc_ins];

% Initial values for signal generator
b_acc = [0.1 0.3 -0.1]'; % IMU biases
b_ars = [0.05 0.1 -0.05]';
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	% Initial states 

% ==============================================================================
% Multirate scheduling
% ==============================================================================
% At most one slow measurement is processed per fast time step.
if f_pos > f_s
    error('The aiding frequency f_pos must satisfy f_pos <= f_s.');
end

slowIndex = 0; % Zero-based index of the next nominal slow measurement

% Tolerance for floating-point comparisons of coincident sample times
tol = 10 * eps(max(1, T_final));

% ==============================================================================
% Time and data initialization
% ==============================================================================
t = 0:h:T_final;
nTimeSteps = length(t);

% Include measurements at t = 0 and, when applicable, t = T_final
N_slow = floor(T_final/h_pos) + 1;

simdata = zeros(nTimeSteps,6); % Pre-allocate table for simulation data
posdata = zeros(N_slow,2);     % Pre-allocate table for position data
positionIndex = 0;

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
for i=1:nTimeSteps
    
    % INS signal generator
    [x, f_imu, ~] = insSignal(x, h, t(i), mu, m_ref);
    phi = x(10); % roll angle
    theta = x(11); % pitch angle

    % Determine whether a new slow pressure measurement is available
    newSlowMeasurement = t(i) + tol >= slowIndex * h_pos;

    if newSlowMeasurement

        slowIndex = slowIndex + 1;
        positionIndex = positionIndex + 1;

        % Pressure aiding
        p = p_0 + rho * g * x(3) + sigma_p * randn;
        z_meas = (p - p_0) / (rho * g);
        posdata(positionIndex,:) = [t(i), z_meas];

        [x_ins, P_prd] = ins_heave( ...
            x_ins, P_prd, h, Qd, Rd, T_acc, mu, rho, ...
            f_imu, phi, theta, p_0, p);

    else

        % No new pressure measurement
        [x_ins, P_prd] = ins_heave( ...
            x_ins, P_prd, h, Qd, Rd, T_acc, mu, rho, ...
            f_imu, phi, theta);
    end

    % Store simulation data in a table (for testing)
    simdata(i,:) = [x(3) x(6) x(9) x_ins'];

end

% Remove unused preallocated rows
posdata = posdata(1:positionIndex,:);

% ==============================================================================
%% PLOTS  
% ==============================================================================
x     = simdata(:,1:3); % High-rate IMU data
x_hat = simdata(:,4:6); 

t_m = posdata(:,1);     % Slow-rate measurements
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

% ==============================================================================
%% DISPLAY DATA
% ==============================================================================
    function displayMethod(cityName, f_s, f_pos);
        disp('-------------------------------------------------------------------');
        disp('MSS toolbox: Error-state Kalman filter (ESKF) for heave estimation');
        disp(['INS aided by pressure measurements at ',num2str(f_pos), ' Hz']);
        disp(['IMU measurements (specific force) at ',num2str(f_s),' Hz']);
        disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
        disp('-------------------------------------------------------------------');
        disp('Simulating...');
    end

end