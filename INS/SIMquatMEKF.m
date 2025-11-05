% SIMquatMEKF is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates the Multiplicative Extended Kalman Filter (MEKF) for 
% quaternion attitude estimation (Markley and Crassidis, 2014). The MEKF uses
% high-rate inertial measurements from a 9-DOF inertial measurement unit 
% (IMU). The MEKF can be called either as a corrector (with new measurements) 
% or as a predictor (without new measurements). The corrector runs at the 
% slowest IMU sensor rate (typically, f_slow is 50 to 100 Hz for a commercial IMU 
% magnetometer), while the predictor runs at the high-rate IMU measurement 
% frequency (typically, f_fast is 500 to 2000 Hz). It is also possible to use a 
% scalar compass measurement (megnetic compass, gyrocompass, GNSS compass, etc.)
% instead of the 3-axis magnetometer measurements.
%
% Dependencies:
%   quatMEKF.m  - MEKF using reference vectors
%   magneticField.m - Magnetic field vectors for different cities%
% 
% See also: SIMquatObserver.m (nonlinear quaternion.based obsever)
%  
% References:
%   Crassidis, J. L., F. L. Markley and Y. Cheng (2007). Survey of Nonlinear 
%       Attitude Estimation Methods. Journal of Guidance, Control and Dynamics 
%       30(1), 12-28.
%
%   Markley, F. L. and J. L. Crassidis (2014). Fundamentals of Spacecraft 
%       Attitude Determination and Control, Volume 33 of Space Technology 
%       Library. Springer-Verlag, New York.
%
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 2nd Edition, Wiley.
%
% Author: Thor I. Fossen
% Date: 2025-11-05
% Revisions: 

clearvars;

% ==============================================================================
% Simulation parameters
% ==============================================================================
T_final = 150; % Final simulation time (s)
f_fast = 1000; % High-rate IMU measurement frequency (Hz)
f_slow = 100; % Low-rate magnetometer/compass measurement frequency (Hz)

% Sampling times in seconds
h_fast = 1/f_fast; % State propagation	 
h_slow = 1/f_slow; % Corrector

% ==============================================================================
% Observer initialization
% ==============================================================================
headingFlag = 1; % 1 for 3-axis magnetometer, 2 for scalar compass

Qd = 1 * diag([1 1 1 1 1 1]); % Gibbs vector and ARS bias noise covariance matrix
Rd = 1 * eye(6); % Gravity and magnetic field/compass measurement covariance matrix
P_prd = 1 * eye(6);

% Initialization of state vectors and covariance matrix
quat_prd = [1 0 0 0]'; 
b_ars_prd = [0 0 0]';

% ==============================================================================
% Initialization of INS signal generator
% ==============================================================================
% Magnetic field and latitude for city #1, see magneticField.m
[m_ref, l, mu, cityName] = magneticField(4);

% IMU biases
b_ars = [0.05 0.05 -0.1]';
   
% Initial values for signal generator
x = [zeros(1,6) zeros(1,3) zeros(1,3) b_ars']';	        

% ==============================================================================
% Display simulation parameters
% ==============================================================================
disp('-------------------------------------------------------------------');
disp('MSS toolbox: Multiplicative Extended Kalman Filter (MEKF) for attitude estimation');
disp(['IMU inertial measurements (specific force and ARS) at ',num2str(f_fast),' Hz']);
disp(['IMU magnetic field/compass measurements at ',num2str(f_slow),' Hz']);
disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
disp('-------------------------------------------------------------------');
disp('Simulating...');

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
t = 0:h_fast:T_final;                % Fast time vector
k = 0;                               % Slow-sample index
next_meas_time = 0;                  % Time for next corrector call
N_slow = floor(T_final/h_slow) + 1;  % Number of slow samples

% Pre-allocate slow-rate storage (1 row per slow sample)
% Columns: [phi theta psi  b_ars(3)  quat_prd(4)  b_ars_prd(3)]
simdata = zeros(N_slow, 13); 

for i=1:length(t)
    
    % INS signal generator, using test signal no. 2 
    [x, f_imu, w_imu, m_imu] = insSignal(x, h_fast, t(i), mu, m_ref, 2);
    phi = x(10); 
    theta = x(11);
    psi = x(12);
    b_ars = x(13:15);

    % Observer correction step runs at slow time
    if t(i) + 1e-10 >= next_meas_time
        switch headingFlag
            case 1 % Heading from magnetometer measurement
                imu_meas = [f_imu' w_imu' m_imu'];
            case 2 % Heading from compass measurement
                imu_meas = [f_imu' w_imu' psi];
        end
        
        next_meas_time = next_meas_time + h_slow;
        do_correct = true;

    else % No magnetometer/compass measurement this fast step
        imu_meas  = [f_imu' w_imu'];
        do_correct = false;
    end

    [quat_prd, b_ars_prd, P_prd] = quatMEKF( ...
        quat_prd, b_ars_prd, P_prd, h_fast, Qd, Rd, m_ref, imu_meas);
       
    % Store slow-rate data 
    if do_correct
        k = k + 1;
        t_meas(k,1) = t(i);
        simdata(k,:) = [phi theta psi  b_ars'  quat_prd'  b_ars_prd'];
    end

end

% ==============================================================================
%% PLOTS         
% ==============================================================================
t_slow = (0:h_slow:T_final)';
phi = simdata(:,1); 
theta = simdata(:,2); 
psi = simdata(:,3); 
b_ars = simdata(:,4:6);  
quat_prd = simdata(:,7:10); 
b_ars_prd = simdata(:,11:13); 

phi_prd = zeros(N_slow,1); 
theta_prd = zeros(N_slow,1); 
psi_prd = zeros(N_slow,1);
for i = 1:N_slow
    [phi_prd(i), theta_prd(i), psi_prd(i)] = q2euler(quat_prd(i,:));
end

% Figure 1
figure(1); figure(gcf);
subplot(311)
plot(t_slow, rad2deg(phi),t_slow, rad2deg(phi_prd)); 
xlabel('Time (s)'),title('Roll angle [deg]'),grid
legend('\phi', '\phi (estimate)');

subplot(312)
plot(t_slow, rad2deg(theta),t_slow, rad2deg(theta_prd)); 
xlabel('Time (s)'),title('Pitch angle [deg]'),grid
legend('\theta', '\theta (estimate)');

subplot(313)
plot(t_slow, rad2deg(psi),t_slow, rad2deg(psi_prd)); 
xlabel('Time (s)'),title('Yaw angle [deg]'),grid
legend('\psi', '\psi (estimate)');

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

% Figure 2
figure(2); figure(gcf);
subplot(311)
plot(t_slow, b_ars(:,1),t_slow, b_ars_prd(:,1)); 
xlabel('Time (s)'),title('ARS roll bias'),grid
legend('b_{x,ars}', 'b_{x,ars} (estimate)');

subplot(312)
plot(t_slow, b_ars(:,2),t_slow, b_ars_prd(:,2)); 
xlabel('Time (s)'),title('ARS pitch bias'),grid
legend('b_{y,ars}', 'b_{y,ars} (estimate)');

subplot(313)
plot(t_slow, b_ars(:,3),t_slow, b_ars_prd(:,3)); 
xlabel('Time (s)'),title('ARS yaw bias'),grid
legend('b_{z,ars}', 'b_{z,ars} (estimate)');

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

