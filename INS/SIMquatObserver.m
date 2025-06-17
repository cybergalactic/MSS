% SIMquatObserver is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates the nonlinear quaternion-based attitude observer of  
% Mahony et al. (2008) using the vector-cross product representation of Grip 
% et al. (2013) with attitude-rate sensor (ARS) bias estimation. The observer 
% uses high-rate inertial measurements from a 9-DOF inertial measurement unit 
% (IMU). The observer can be called either as a corrector (with new 
% measurements) or as a predictor (without new measurements). Additionally, 
% the magnetometer can operate at a slower rate (typically 100 Hz) compared to 
% the high-rate specific force and ARS measurements (typically 1000 Hz). It 
% is also possible to use a scalar compass measurement.
%
% Dependencies:
%   quatObserver.m  - Nonlinear attiitude observer using reference vectors
%   magneticField.m - Magnetic field vectors for different cities
%  
% References:
%   H. F. Grip, T. I. Fossen, T. A. Johansen, and A. Saberi (2013). 
%       Nonlinear Observer for GNSS-Aided Inertial Navigation with 
%       Quaternion-Based Attitude Estimation. American Control Conference, 
%       Washington DC, USA, IEEE Xplore, pp. 272-279. 
%       doi.org/10.1109/ACC.2013.6579849
%
%  R. Mahony, T. Hamel and J.-M. Pflimlin (2008). Nonlinear Complementary 
%       Filters on the Special Orthogonal Group. IEEE Trans. on Aut. Control 53(5)

%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 2nd Edition, Wiley.
%
% Author: Thor I. Fossen
% Date: 2024-08-20
% Revisions: 
%   2025-06-17 Modified to accept compass measurements.

% ==============================================================================
% Simulation parameters
% ==============================================================================
T_final = 200; % Final simulation time (s)
f_fast = 1000; % High-rate IMU measurement frequency (Hz)
f_slow = 100; % Low-rate magnetometer/compass measurement frequency (Hz)

% Sampling times in seconds
h_fast = 1/f_fast; % State propagation	 
h_slow = 1/f_slow; % Corrector

% ==============================================================================
% Observer initialization
% ==============================================================================
headingFlag = 1; % 1 for 3-axis magnetometer, 2 for scalar compass
coningSculling = 0; % 0 for no compensation, 1 for compensation of coning and sculling

switch headingFlag
    case 1
        k1 = 1000; % Gain for specific force measurement vector
        k2 = 500; % Gain for magnetic field measurement vector
        Ki = 0.2 * diag([ 1 1 1 ]); % Integral gain matrix for ARS bias estimation
    case 2
        k1 = 1500; % Gain for specific force measurement vector
        k2 = 50; % Gain for compass measurement 
        Ki = 0.2 * diag([ 1 1 0.5 ]); % Integral gain matrix for ARS bias estimation
end

% Initialization of observer states
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
disp('MSS toolbox: Nonlinear quaternion-based attitude observer');
disp(['IMU inertial measurements (specific force and ARS) at ',num2str(f_fast),' Hz']);
disp(['IMU magnetic field/compass measurements at ',num2str(f_slow),' Hz']);
disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
if coningSculling == 0
    disp('No coning and sculling compensation');
else
    disp('Coning and sculling compensation by the midtpoint method');
end
disp('-------------------------------------------------------------------');
disp('Simulating...');

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
t = 0:h_fast:T_final; % Time vector from 0 to T_final          
nTimeSteps = length(t); % Number of time steps
simdata = zeros(nTimeSteps,13); % Pre-allocate table for simulation data 

for i=1:nTimeSteps
    
    % INS signal generator, using test signal no. 2 
    [x, f_imu, w_imu, m_imu] = insSignal(x, h_fast, t(i), mu, m_ref, 2);
    phi = x(10); 
    theta = x(11);
    psi = x(12);
    b_ars = x(13:15);

    % Observer correction step runs at slow time
    if mod( t(i), h_slow ) == 0
        switch headingFlag
            case 1 % Heading from magnetometer measurement
                imu_meas = [f_imu' w_imu' m_imu'];
            case 2 % Heading from compass measurement
                imu_meas = [f_imu' w_imu' psi];
        end
    else  % No magnetometer/compass measurement
        imu_meas = [f_imu' w_imu'];
    end

    [quat_prd, b_ars_prd] = quatObserver( ...
        quat_prd, b_ars_prd, h_fast, Ki, k1, k2, m_ref, imu_meas, coningSculling);
       
    % Store simulation data in a table 
    simdata(i,:) = [phi theta psi b_ars' quat_prd' b_ars_prd']; 
    
end

% ==============================================================================
%% PLOTS         
% ==============================================================================
phi = simdata(:,1); 
theta = simdata(:,2); 
psi = simdata(:,3); 
b_ars = simdata(:,4:6);  
quat_prd = simdata(:,7:10); 
b_ars_prd = simdata(:,11:13); 

phi_prd = zeros(length(t),1); 
theta_prd = zeros(length(t),1); 
psi_prd = zeros(length(t),1);
for i = 1:length(t)
    [phi_prd(i), theta_prd(i), psi_prd(i)] = q2euler(quat_prd(i,:));
end

% Figure 1
figure(1); figure(gcf);
subplot(311)
plot(t, rad2deg(phi),t, rad2deg(phi_prd)); 
xlabel('Time (s)'),title('Roll angle [deg]'),grid
legend('\phi', '\phi_{prd}');

subplot(312)
plot(t, rad2deg(theta),t, rad2deg(theta_prd)); 
xlabel('Time (s)'),title('Pitch angle [deg]'),grid
legend('\theta', '\theta_{prd}');

subplot(313)
plot(t, rad2deg(psi),t, rad2deg(psi_prd)); 
xlabel('Time (s)'),title('Yaw angle [deg]'),grid
legend('\psi', '\psi_{prd}');

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

% Figure 2
figure(2); figure(gcf);
subplot(311)
plot(t, b_ars(:,1),t, b_ars_prd(:,1)); 
xlabel('Time (s)'),title('ARS roll bias'),grid
legend('bx_{ars}', 'bx_{ars_prd}');

subplot(312)
plot(t, b_ars(:,2),t, b_ars_prd(:,2)); 
xlabel('Time (s)'),title('ARS pitch bias'),grid
legend('by_{ars}', 'by_{ars_prd}');

subplot(313)
plot(t, b_ars(:,3),t, b_ars_prd(:,3)); 
xlabel('Time (s)'),title('ARS yaw bias'),grid
legend('bz_{ars}', 'bz_{ars_prd}');

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

