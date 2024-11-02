function SIMaidedINSeuler()
% SIMaidedINSeuler is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates an Inertial Navigation System (INS) aided by position 
% measurements using the Error-State Kalman Filter (ESKF). The attitude is 
% parametrized using Euler angles (Fossen, 2021, Chapter 14.4).  
% 
% The position measurement frequency f_pos (typically 5 Hz) can be chosen smaller 
% or equal to the sampling frequency f_s (typically 1000 Hz), which is equal to 
% the Inertial Measurement Unit (IMU) measurement frequency. 
%
% Dependencies:
%   ins_euler.m     - Feedback ESKF for INS aided by position measurements 
%                     y_pos and compass measurements y_psi. The velocity aiding 
%                     signal y_vel is optionally. 
%   ins_ahrs.m      - Feedback ESKF for INS aided by position measurements 
%                     y_pos and attitude measurements (roll, pitch and yaw 
%                      angles). The velocity aiding signal y_vel is optionally.
%   magneticField.m - Magnetic field vectors for different cities.
%  
% References:
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2021-04-26
% Revisions:
%   2024-08-20 : Using the updated insSignal.m generator.
%   2024-11-02 : Improved logic for slow position data

%% USER INPUTS
T_final = 100; % Final simulation time (s)
f_s    = 1000; % Sampling frequency equals IMU measurement frequency (Hz)
f_pos = 5; % Position measurement frequency (Hz)

% Sampling times
h  = 1/f_s; 	 
h_pos = 1/f_pos; 

% Initialization of the INS signal generator
[m_ref, ~,mu,cityName] = magneticField(1); % Magntic field and latitude for city #1
b_acc = [0.1 0.3 -0.1]'; % IMU biases
b_ars = [0.05 0.1 -0.05]';
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']'; % Initial states for signal generator        

% Display simulation options
[attitudeFlag, velFlag] = displayMethod(cityName);

% Initialization of ESKF covariance matrix
P_prd = eye(15);

if (attitudeFlag == 1) % Compass    
    % Process noise weights: vel, acc_bias, w_nb, ars_bias
    Qd = diag([0.1 0.1 0.1  0.001 0.001 0.001  0.1 0.1 0.1  0.001 0.001 0.001]);   
    if (velFlag == 1)
        % Position and compass aiding
        Rd = diag([0.1 0.1 0.1  1 1 1  0.001]); % pos, acc, compass
    else % velFlag == 2
        % Position/velocity aiding + compass
        Rd = diag([1 1 1  1 1 1  1 1 1  0.001]); % pos, vel, acc, psi
    end
else % attitudeFlag == 2 (AHRS)
    if (velFlag == 1) 
       Rd = diag([1 1 1  1 1 1]); % pos, euler_angles
       Qd = diag([1 1 1  1 1 1  10 10 10  0.01 0.01 0.01]);
    else 
       Rd = diag([10 10 10 1 1 1 0.1 0.1 0.1]); % pos, vel, euler_angles
       Qd = diag([1 1 1  1 1 1  0.1 0.1 0.1  0.01 0.01 0.01]); 
    end
end

% Initialization of INS states
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
theta_ins = [0, 0, 0]';
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; theta_ins; b_ars_ins];

% Time vector initialization
t_slow = 0; % Initialize the time for the next slow measurement
t = 0:h:T_final; % Time vector from 0 to T_final          
nTimeSteps = length(t); % Number of time steps

%% MAIN LOOP
simdata = zeros(nTimeSteps,30); % Pre-allocate table for simulation data
posdata = zeros(floor(T_final * f_pos), 4); % Pre-allocate table for position data
pos_index = 0; % Initialize index for posdata

for i=1:nTimeSteps
    
    % INS signal generator
    [x, f_imu, w_imu] = insSignal(x, h, t(i), mu, m_ref);
    y_psi = x(12);
    y_ahrs = x(10:12);
    
    % Positions measurements are slower than the sampling time
    if t(i) > t_slow
        % Aiding
        pos_index = pos_index + 1; 
        y_pos = x(1:3) + 0.05 * randn(3,1); % Position measurements
        y_vel = x(4:6) + 0.01 * randn(3,1); % Optionally velocity meas.
        posdata(pos_index, :) = [t(i), y_pos']; % Store position measurements 

        if (attitudeFlag == 1) % Compass
            if (velFlag == 1)
                % Position aiding + compass aiding
                [x_ins,P_prd] = ins_euler(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos);
            else
                % Position/velocity aiding + compass aiding
                [x_ins,P_prd] = ins_euler(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos,y_vel);
            end
        else % attitudeFlag == 2 (AHRS)
            if (velFlag == 1)
                [x_ins,P_prd] = ins_ahrs(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs,y_pos);
            else
                [x_ins,P_prd] = ins_ahrs(...
                    x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs,y_pos,y_vel);
            end
        end
        
        % Update the time for the next slow position measurement
        t_slow = t_slow + h_pos; 

    else  
        % No aiding
        if (attitudeFlag == 1) 
            % Compass
            [x_ins,P_prd] = ins_euler(x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi);
        else 
            % AHRS
            [x_ins,P_prd] = ins_ahrs(x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_ahrs);
        end
    end

    % Store simulation data in a table (for testing)
    simdata(i,:) = [x' x_ins'];

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Get screen dimensions
legendSize = 12;
         
x     = simdata(:,1:15); % High-rate IMU data
x_hat = simdata(:,16:30); 

t_m = posdata(:,1); % Slow-rate position measurements
y_m = posdata(:,2:4);

figure(1); 
if ~isoctave;set(gcf,'Position',[1,1,0.4*scrSz(3),scrSz(4)]);end

subplot(311)
h1 = plot(t_m,y_m,'xr'); hold on;
h2 = plot(t,x_hat(:,1:3),'b'); hold off;
xlabel('time (s)'),title('Position [m]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_pos), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

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
set(findall(gcf,'type','legend'),'FontSize',legendSize)

figure(2); 
if ~isoctave;set(gcf,'Position',[0.4*scrSz(3),1,0.4*scrSz(3),scrSz(4)]);end

subplot(211)
h1 = plot(t,rad2deg(x(:,10:12)),'r'); hold on;
h2 = plot(t,rad2deg(x_hat(:,10:12)),'b'); hold off;
xlabel('time (s)'),title('Angle [deg]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(212)
h1 = plot(t,x(:,13:15),'r'); hold on;
h2 = plot(t,x_hat(:,13:15),'b'); hold off;
xlabel('time (s)'),title('ARS bias'),grid
legend([h1(1),h2(1)],['True ARS bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% RADIO BUTTONS, FLAGS AND DISPLAY
function [attitudeFlag, velFlag] = displayMethod(cityName)

    f = figure('Position', [400, 400, 400, 300], 'Name', 'Strapdown Aided INS', 'MenuBar', 'none', 'NumberTitle', 'off', 'WindowStyle', 'modal');

    % Add button group for control methods
    bg1 = uibuttongroup('Parent', f, 'Position', [0.02 0.65 0.96 0.3], 'Title', 'Attitude Aiding','FontSize',14,'FontWeight','bold');
    radio1 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',13, 'String', 'Compass', 'Position', [10 40 500 30], 'Tag', '1');
    radio2 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',13, 'String', 'Attitude and heading reference system (AHRS)', 'Position', [10 10 500 30], 'Tag', '2');
    set(radio1, 'Value', 1); % Set default value

    % Add button group for velocity aiding options
    bg2 = uibuttongroup('Parent', f, 'Position', [0.02 0.35 0.96 0.3], 'Title', 'Velocity Aiding','FontSize',14,'FontWeight','bold');
    radio3 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 13, 'String', 'No Velocity Aiding', 'Position', [10 35 500 30], 'Tag', '1');
    radio4 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 13, 'String', 'Velocity Aiding', 'Position', [10 5 500 30], 'Tag', '2');
    set(radio3, 'Value', 1); % Set default value

    % Add OK button to confirm selections
    uicontrol('Style', 'pushbutton', 'String', 'OK', 'FontSize', 13, 'Position', [20 30 100 40], 'Callback', @(src, evt) uiresume(f));

    uiwait(f); % wait for uiresume to be called on figure handle

    % Determine which attitude method was selected
    if get(radio1, 'Value') == 1
        attitudeFlag = str2double(get(radio1, 'Tag'));
    else
        attitudeFlag = str2double(get(radio2, 'Tag'));
    end

    % Determine if velocity aiding was selected
    if get(radio3, 'Value') == 1
        velFlag  = str2double(get(radio3, 'Tag'));
    else
        velFlag  = str2double(get(radio4, 'Tag'));
    end

    close(f);  % close the figure after obtaining the selections

    disp('-------------------------------------------------------------------');
    disp('MSS toolbox: Error-state (indirect) feedback Kalman filter');
    disp('Attitude parametrization: Euler angles');
    if (velFlag == 1)
        disp(['INS aided by position at ',num2str(f_pos), ' Hz']);
    else
        disp(['INS aided by position and velocity at ',num2str(f_pos),' Hz']);
    end
    disp(['IMU measurements (specific force and ARS) at ',num2str(f_s),' Hz']);
    if (attitudeFlag == 1)
        disp(['Compass measurements at ',num2str(f_s), ' Hz']);
    else
        disp(['Three-axis AHRS measurements at ',num2str(f_s), ' Hz']);
    end
    disp(['Magnetic field reference vector for ', cityName, ' (>> type magneticField)']);
    disp('-------------------------------------------------------------------');
    disp('Simulating...');

end

end

