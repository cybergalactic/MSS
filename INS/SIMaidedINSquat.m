function SIMaidedINSquat()
% SIMaidedINSquat is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates an Inertial Navigation System (INS) aided by position 
% measurements using the Error-State Kalman Filter (ESKF). The attitude is 
% parametrized using unit quaternions and the error states are represented 
% by Gibbs vector in a Multiplicative Extended Kalman Filter (MEKF) 
% (Fossen, 2021, Chapter 14.4).  
%
% The ESKF uses high-rate inertial measurements from a 9-DOF inertial measurement
% unit (IMU). The ESKF can be called either as a corrector (with new measurements) 
% or as a predictor (without new measurements). The high-rate IMU frequency is
% equal to the sampling frequency f_fast (typically 1000 Hz). Additionally, the 
% magnetometer or compass can operate at a slower rate (typically 100 Hz). The 
% position measurement frequency f_slow (typically 5 Hz) must be smaller or equal 
% to the sampling frequency f_fast. 
%
% Dependencies:
%   ins_mekf_psi.m  - Feedback ESKF for INS aided by position measurements 
%                     y_pos and compass measurements y_psi. The velocity 
%                     aiding signal y_vel is optionally. 
%   ins_mekf.m      - Feedback ESKF for INS aided by position measurements 
%                     y_pos and magnetometer measurements y_mag. The velocity 
%                     aiding signal y_vel is optionally.
%   magneticField.m - Magnetic field vectors for different cities.
%  
% References:
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2024-04-26
% Revisions:
%   2024-08-20 : Using the updated insSignal.m generator
%   2024-09-09 : New logic for slow position and magnetometer/compass measurements
%   2025-11-11 : Added application modes and pseudomeasurement for average
%                sea level
clearvars;

% ==============================================================================
% Simulation parameters
% ==============================================================================
T_final = 100; % Final simulation time (s)
f_fast = 1000; % High-rate IMU meaurement frequency (Hz)
f_mag = 100; % Magnetometer measurement frequency (Hz)
f_slow = 5; % Slow-rate position measurement frequency (Hz)

% Sampling times in seconds
h  = 1/f_fast; 	
h_slow = 1/f_slow; 
h_mag = 1/f_mag;

% ==============================================================================
% Initialization of the INS signal generator
% ==============================================================================
[m_ref, ~, mu, cityName] = magneticField(1); % Magntic field and latitude for city #1
b_acc = [0.1 0.3 -0.1]'; %  IMU accelerometer bias
b_ars = [0.05 0.1 -0.05]'; %  % IMU ARS bias
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	% Initial states 

% Display simulation options
[attitudeFlag,velFlag,applicationFlag,pseudoFlag] = displayMethod(cityName); 

% ==============================================================================
% Initialization of ESKF covariance matrices
% ==============================================================================
P_prd = eye(15);

% Process noise weights: vel, acc_bias, w, ars_bias
Qd = diag([0.01 0.01 0.01  0.01 0.01 0.01  0.1 0.1 0.1  0.001 0.001 0.001]);
if pseudoFlag % Additional sea level state
    Qd = blkdiag(Qd, 0.01);
    P_prd = blkdiag(P_prd, 1);
end

% Covariance weights for position, velocity, magnetometer and compass measurements
Rd.position = 1 * diag([1 1 1]);
Rd.velocity = 0.1 * diag([1 1 1]);
Rd.magnetometer = 0.01 * diag([1 1 1]);
Rd.compass = 0.01;

% Application-based aiding
Rd.gravityRefVector = 1 * diag([1 1 1]); % Gravity reference vector (hovering/stationkeeping)
Rd.pseduoMeasSeaLevel = 0.01; % Pseduomeasurement for average sea level

% Configuration flags used by ins_mekf.m
Rd.applicationFlag = applicationFlag;
Rd.pseudoFlag = pseudoFlag;

% ==============================================================================
% Initialization of the INS states
% ==============================================================================
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
q_ins = euler2q(0, 0, 0);
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];
if pseudoFlag, x_ins = [x_ins; 0]; end % Additional sea level state

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps
simdata = zeros(nTimeSteps,31); % Pre-allocate table for simulation data
if pseudoFlag, simdata = zeros(nTimeSteps,32); end
posdata = zeros(floor(T_final * f_slow), 4); % Pre-allocate table for position data
pos_index = 0; % Initialize index for posdata

for i=1:nTimeSteps
    
    % INS signal generator 
    [x, f_imu, w_imu, m_imu] = insSignal(x, h, t(i), mu, m_ref);
   
     % IMU magnetometer and compass measurements are slower than the sampling time
    if abs( mod(t(i), h_mag) ) < 1e-10
        imu_meas = [f_imu' w_imu' m_imu']; % 9-DOF IMU measurements
        y_psi = x(12); % Compass measurement
    else 
        imu_meas = [f_imu' w_imu']; % No magnetometer measurements
        y_psi = []; % No compass measurement
    end

    % Position measurements are slower than the sampling time
    if abs(mod(t(i), h_slow)) < 1e-10
        % AIDING
        pos_index = pos_index + 1;
        y_pos = x(1:3) + 0.05 * randn(3,1);
        y_vel = x(4:6) + 0.01 * randn(3,1);
        posdata(pos_index, :) = [t(i), y_pos'];

        if attitudeFlag == 2   
            % MAGNETOMETER 
            if ~velFlag
                % Position aiding + magnetometer
                [x_ins,P_prd] = ins_mekf(...
                    x_ins,P_prd,mu,h,Qd,Rd,imu_meas,m_ref,y_pos);
            else
                % Position/velocity aiding + magnetometer
                [x_ins,P_prd] = ins_mekf(...
                    x_ins,P_prd,mu,h,Qd,Rd,imu_meas,m_ref,y_pos,y_vel);
            end

        else  
            % COMPASS
            if ~velFlag
                % Position aiding + compass
                [x_ins,P_prd] = ins_mekf_psi(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos);
            else
                % Position/velocity aiding + compass
                [x_ins,P_prd] = ins_mekf_psi(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos,y_vel);
            end
        end

    else
        % NO AIDING 
        if attitudeFlag == 2
            % Magnetometer
            [x_ins,P_prd] = ins_mekf(...
                x_ins,P_prd,mu,h,Qd,Rd,[f_imu', w_imu', m_imu'],m_ref);
        else
            % Compass
            [x_ins,P_prd] = ins_mekf_psi(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi);
        end
    end

    % Store simulation data in a table (for testing)
    simdata(i,:) = [x' x_ins']; 
    
end

% ==============================================================================
% PLOTS
% ==============================================================================
scrSz = get(0, 'ScreenSize'); % Get screen dimensions
legendSize = 10;
colors = {'b','g','k'};

x     = simdata(:,1:15); 
x_hat = simdata(:,16:31); 

phi = zeros(length(t),1);theta = zeros(length(t),1);psi = zeros(length(t),1);
for i = 1:length(t)
    [phi(i), theta(i), psi(i)] = q2euler(x_hat(i,10:13));
end
Theta = [phi theta psi];

t_m = posdata(:,1); % Slow-rate position measurements
y_m = posdata(:,2:4);

% Figure 1
figure(1); 
if ~isoctave; set(gcf,'Position',[1,1,0.4*scrSz(3),scrSz(4)]); end

subplot(311)
hTrue = plot(t_m, y_m, 'xr');
hold on
hX = plot(t, x_hat(:,1), colors{1});
hY = plot(t, x_hat(:,2), colors{2});
hZ = plot(t, x_hat(:,3), colors{3});
hAll = [hX, hY, hZ, hTrue(1)];
hold off
xlabel('Time [s]'); title('Position [m]'); grid on
labels = { ...
    ['Estimate x_N at ', num2str(f_fast), ' Hz'], ...
    ['Estimate y_E at ', num2str(f_fast), ' Hz'], ...
    ['Estimate z_D at ', num2str(f_fast), ' Hz'], ...
    ['Position measurements at ', num2str(f_slow), ' Hz'] };
legend(hAll, labels);

subplot(312)
hTrue = plot(t,x(:,4:6),'r');
hold on
hX = plot(t, x_hat(:,4), colors{1});
hY = plot(t, x_hat(:,5), colors{2});
hZ = plot(t, x_hat(:,6), colors{3});
hAll = [hX, hY, hZ, hTrue(1)];
hold off
xlabel('Time [s]'),title('Velocity [m/s]'),grid
labels = { ... 
    ['Estimate v_N at ', num2str(f_fast), ' Hz'], ...
    ['Estimate v_E at ', num2str(f_fast), ' Hz'], ...
    ['Estimate v_D at ', num2str(f_fast), ' Hz'], ...
    ['True velocity at ', num2str(f_fast), ' Hz']};
legend(hAll, labels);

subplot(313)
hTrue = plot(t,x(:,7:9),'r');
hold on
hX = plot(t, x_hat(:,7), colors{1});
hY = plot(t, x_hat(:,8), colors{2});
hZ = plot(t, x_hat(:,9), colors{3});
hAll = [hX, hY, hZ, hTrue(1)];
hold off
xlabel('Time [s]'),title('Acceleration bias [m/s^2]'),grid
labels = { ... 
    ['Estimate b_{x, acc} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate b_{y, acc} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate b_{z, acc} at ', num2str(f_fast), ' Hz'], ...
    ['True acceleration bias at ', num2str(f_fast), ' Hz']};
legend(hAll, labels);

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

% Figure 2
figure(2); 
if ~isoctave;set(gcf,'Position',[0.4*scrSz(3),1,0.4*scrSz(3),scrSz(4)]); end

subplot(211)
hTrue = plot(t,rad2deg(x(:,10:12)),'r');
hold on
hX = plot(t,rad2deg(Theta(:,1)), colors{1});
hY = plot(t,rad2deg(Theta(:,2)), colors{2});
hZ = plot(t,rad2deg(Theta(:,3)), colors{3});
hAll = [hX, hY, hZ, hTrue(1)];
hold off
xlabel('Time [s]'),title('Euler angles [deg]'),grid
labels = { ... 
    ['Estimate \phi at ', num2str(f_fast), ' Hz'], ...
    ['Estimate \theta at ', num2str(f_fast), ' Hz'], ...
    ['Estimate \psi at ', num2str(f_fast), ' Hz'], ...
    ['True Euler angles at ', num2str(f_fast), ' Hz']};
legend(hAll, labels);

subplot(212)
hTrue = plot(t,rad2deg(x(:,13:15)),'r');
hold on
hX = plot(t,rad2deg(x_hat(:,14)), colors{1});
hY = plot(t,rad2deg(x_hat(:,15)), colors{2});
hZ = plot(t,rad2deg(x_hat(:,16)), colors{3});
hAll = [hX, hY, hZ, hTrue(1)];
hold off
xlabel('Time [s]'),title('Angular rate bias [deg/s]'),grid
labels = { ... 
    ['Estimate a_{x, ars} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate a_{y, ars} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate a_{z, ars} at ', num2str(f_fast), ' Hz'], ...
    ['True ARS bias at ', num2str(f_fast), ' Hz']};
legend(hAll, labels);

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

% ==============================================================================
%% RADIO BUTTONS, FLAGS AND DISPLAY
% ==============================================================================
function [attitudeFlag,velFlag,applicationFlag,pseudoFlag] = displayMethod(cityName)

    f = figure('Position', [400, 300, 550, 600], 'Name', 'Strapdown Aided INS', 'MenuBar', 'none', 'NumberTitle', 'off', 'WindowStyle', 'modal');

    % Add button group for control methods
    bg1 = uibuttongroup('Parent', f, 'Position', [0.02 0.67 0.96 0.3], 'Title', 'Compass Aiding','FontSize',14,'FontWeight','bold');
    radio1 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',12, 'String', 'Compass', 'Position', [10 120 500 30], 'Tag', '1');
    radio2 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',12, 'String', 'Magnetometer', 'Position', [10 90 500 30], 'Tag', '2');
    set(radio2, 'Value', 1); % Set default value

    % Add button group for velocity aiding options
    bg2 = uibuttongroup('Parent', f, 'Position', [0.02 0.49 0.96 0.3], 'Title', 'Velocity Aiding','FontSize',14,'FontWeight','bold');
    radio3 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'No Velocity Aiding', 'Position', [10 120 500 30], 'Tag', 'false');
    radio4 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'Velocity Aiding', 'Position', [10 90 500 30], 'Tag', 'true');
    set(radio3, 'Value', 1);     % Default = NO velcoity measurement

    % Add button group for application-based aiding
    bg3 = uibuttongroup('Parent', f, 'Position', [0.02 0.31 0.96 0.3], 'Title', 'Application-Based Aiding Mode','FontSize',14,'FontWeight','bold');
    radio5 = uicontrol(bg3, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'No gravity Reference Vector Measurement (Large Roll/Pitch Angles)', 'Position', [10 120 500 30], 'Tag', 'false');
    radio6 = uicontrol(bg3, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'Gravity Reference Vector Measurement (Small Roll/Pitch Angles)', 'Position', [10 90 500 30], 'Tag', 'true');
    set(radio5, 'Value', 1);   % Default = NO aiding

    % Add button group for pseudomeasurement constraint
    bg4 = uibuttongroup('Parent', f, 'Position', [0.02 0.12 0.96 0.3], 'Title', 'Pseudomeasurement (Average Sea Surface zⁿ ≈ 0)','FontSize',14,'FontWeight','bold');
    radio7 = uicontrol(bg4, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'No Pseudomeasurement', 'Position', [10 120 500 30], 'Tag', 'false');
    radio8 = uicontrol(bg4, 'Style', 'radiobutton', 'FontSize', 12, 'String', 'Enable Pseudomeasurement', 'Position', [10 90 500 30], 'Tag', 'true');
    set(radio7, 'Value', 1);   % Default = NO pseudomeasurement

    % Add OK button to confirm selections
    uicontrol('Style', 'pushbutton', 'String', 'OK', 'FontSize', 12, 'Position', [20 30 100 40], 'Callback', @(src, evt) uiresume(f));

    uiwait(f); % wait for uiresume to be called on figure handle

    attitudeFlag    = str2double(get(findobj(bg1,'Value',1),'Tag'));
    velFlag         = strcmp(get(findobj(bg2,'Value',1),'Tag'), 'true');
    applicationFlag = strcmp(get(findobj(bg3,'Value',1),'Tag'), 'true');
    pseudoFlag      = strcmp(get(findobj(bg4,'Value',1),'Tag'), 'true');

    close(f);  % close the figure after obtaining the selections

    disp('-------------------------------------------------------------------');
    disp('MSS toolbox: Error-state (indirect) feedback Kalman filter');
    disp('Unit quaternion attitude parametrization (MEKF): 2 x Gibbs vector');
   
    if ~velFlag
        disp(['INS aided by position at ',num2str(f_slow), ' Hz']);
    else
        disp(['INS aided by position and velocity at ',num2str(f_slow),' Hz']);
    end

    disp(['IMU inertial measurements (specific force and ARS) at ',num2str(f_fast),' Hz']);
    
    if (attitudeFlag == 2)
        disp(['IMU magnetometer measurements at ',num2str(f_mag), ' Hz']);
        disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
    else
        disp(['Compass measurements at ',num2str(f_mag), ' Hz']);
    end

    if applicationFlag
        disp('Gravity reference vector (hover/stationkeeping)');
    end  

    if pseudoFlag
        disp('Pseudomeasurement for average sea surface (zⁿ ≈ 0)');
    end 

    disp('-------------------------------------------------------------------');
    disp('Simulating...');

end

end