function SIMaidedINSquat()
% SIMaidedINSquat is compatible with MATLAB and GNU Octave (www.octave.org).
% This script simulates an Inertial Navigation System (INS) aided by position 
% measurements using the Error-State Kalman Filter (ESKF). The attitude is 
% parametrized using unit quaternions and the error states are represented 
% by Gibbs vector in a Multiplicative Extended Kalman Filter (MEKF) 
% (Fossen, 2021, Chapter 14.4).  
% 
% The position measurement frequency f_pos can be chosen smaller 
% or equal to the sampling frequency f_s, which is equal to the Inertial 
% Measurement Unit (IMU) measurement frequency. The ratio between the 
% frequencies must be an integer Z such that:
%
%     Integer:       Z = f_s/f_pos >= 1 
%
% Dependencies:
%   ins_mekf_psi.m  - Feedback ESKF for INS aided by position measurements 
%                     y_pos and compass measurements y_psi. The velocity 
%                     aiding signal y_vel is optionally. 
%   ins_mekf.m      - Feedback ESKF for INS aided by position measurements 
%                     y_pos and magnetometer measurements y_mag. The velocity 
%                     aiding signal y_vel is optionally.
%  
% References:
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2024-04-26
% Revisions:

%% USER INPUTS
N  = 20000;		  % No. of iterations
f_s    = 100;     % Sampling frequency [Hz]
f_pos = 1;        % Position measurement frequency [Hz]

[attitudeFlag, velFlag] = displayMethod();

% Sampling times
h  = 1/f_s; 	 
h_pos = 1/f_pos;  

% IMU biases
b_acc = [0.1 0.3 -0.1]';
b_ars = [0.05 0.1 -0.05]';
   
% Initial values for signal generator
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	        

% Initialization of ESKF covariance matrix
P_prd = eye(15);

% Process noise weights: v, acc_bias, w, ars_bias
Qd = diag([0.01 0.01 0.01 0.01 0.01 0.01 0.1 0.1 0.1 0.001 0.001 0.001]);
   
if (velFlag == 1 && attitudeFlag == 2)
    % Position aiding + magnetometer
    Rd = diag([1 1 1  1 1 1 0.01 0.01 0.01]);  % p, acc, mag
elseif  (velFlag == 2 && attitudeFlag == 2)
    % Position/velocity aiding + magnetometer
    Rd = diag([1 1 1  0.1 0.1 0.1  1 1 1  0.01 0.01 0.01]);  % p, v, acc, mag
elseif (velFlag == 1 && attitudeFlag == 1)
    % Position aiding + compass
    Rd = diag([1 1 1  1 1 1  0.01]);  % p, acc, psi
else
    % Position/velocity aiding + compass
    Rd = diag([1 1 1  1 1 1  1 1 1  0.01]);  % p, vel, acc, psi
end

% Initialization of INS states
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
q_ins = euler2q(0, 0, 0);
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];

% WGS-84 gravity model
mu = 63.4305 * pi / 180;    % lattitude  
g = gravity(mu);  

%% MAIN LOOP
simdata = zeros(N+1,32);                  % Table of simulation data
ydata = [0 x(1:3)'];                      % Table of position measurements

for i=1:N+1
    
    % INS signal generator
    t = (i-1) * h;                        % Time (s)   
    [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t);
    y_psi = x(12);
    
    % Position measurements are slower than the sampling time
    if mod( t, h_pos ) == 0
        
        y_pos = x(1:3) + 0.05 * randn(3,1);   % Position measurements
        y_vel = x(4:6) + 0.01 * randn(3,1);   % Optionally velocity meas.
        ydata = [ydata; t, y_pos'];           % Store position measurements                  
              
        if (velFlag == 1 && attitudeFlag == 2)
            % Position aiding + magnetometer
            [x_ins,P_prd] = ins_mekf(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref,y_pos);

        elseif (velFlag == 2 && attitudeFlag == 2)
            % Position/velocity aiding + magnetometer
            [x_ins,P_prd] = ins_mekf(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref,y_pos,y_vel);

        elseif (velFlag == 1 && attitudeFlag == 1)
            % Position aiding + compass
            [x_ins,P_prd] = ins_mekf_psi(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos);

        else
            % Position/velocity aiding + compass
            [x_ins,P_prd] = ins_mekf_psi(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos,y_vel);
        end
        
    else  % No aiding
        
        if (attitudeFlag == 2)
            % Magnetometer
            [x_ins,P_prd] = ins_mekf(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref);

        else
            % Compass
            [x_ins,P_prd] = ins_mekf_psi(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi);

        end
        
    end
       
    % Store simulation data in a table (for testing)
    simdata(i,:) = [t x' x_ins']; 
    
end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Get screen dimensions
legendSize = 12;

t     = simdata(:,1);           
x     = simdata(:,2:16); 
x_hat = simdata(:,17:32); 

phi = zeros(length(t),1);theta = zeros(length(t),1);psi = zeros(length(t),1);
for i = 1:length(t)
    [phi(i), theta(i), psi(i)] = q2euler(x_hat(i,10:13));
end
Theta = [phi theta psi];

t_m = ydata(:,1);              % Slow position measurements
y_m = ydata(:,2:4);

figure(1); 
if ~isoctave; set(gcf,'Position',[1,1,0.4*scrSz(3),scrSz(4)]); end

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
h1 = plot(t,rad2deg( x(:,10:12) ),'r'); hold on;
h2 = plot(t,rad2deg( Theta ),'b'); hold off;
xlabel('time (s)'),title('Angle [deg]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(212)
h1 = plot(t,x(:,13:15),'r'); hold on;
h2 = plot(t,x_hat(:,14:16),'b'); hold off;
xlabel('time (s)'),title('ARS bias'),grid
legend([h1(1),h2(1)],['True ARS bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',legendSize)

%% RADIO BUTTONS, FLAGS AND DISPLAY
function [attitudeFlag, velFlag] = displayMethod()

    f = figure('Position', [400, 400, 400, 300], 'Name', 'Strapdown Aided INS', 'MenuBar', 'none', 'NumberTitle', 'off', 'WindowStyle', 'modal');

    % Add button group for control methods
    bg1 = uibuttongroup('Parent', f, 'Position', [0.02 0.65 0.96 0.3], 'Title', 'Compass Aiding','FontSize',14,'FontWeight','bold');
    radio1 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',13, 'String', 'Compass', 'Position', [10 40 500 30], 'Tag', '1');
    radio2 = uicontrol(bg1, 'Style', 'radiobutton', 'FontSize',13, 'String', 'Magnetometer', 'Position', [10 10 500 30], 'Tag', '2');

    % Add button group for velocity aiding options
    bg2 = uibuttongroup('Parent', f, 'Position', [0.02 0.35 0.96 0.3], 'Title', 'Velocity Aiding','FontSize',14,'FontWeight','bold');
    radio3 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 13, 'String', 'No Velocity Aiding', 'Position', [10 35 500 30], 'Tag', '1');
    radio4 = uicontrol(bg2, 'Style', 'radiobutton', 'FontSize', 13, 'String', 'Velocity Aiding', 'Position', [10 5 500 30], 'Tag', '2');

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
    disp('Unit quaternion attitude parametrization (MEKF): 2 x Gibbs vector');
    if (velFlag == 1)
        disp(['INS aided by position at ',num2str(f_pos), ' Hz']);
    else
        disp(['INS aided by position and velocity at ',num2str(f_pos),' Hz']);
    end
    disp(['IMU measurements (specific force and ARS) at ',num2str(f_s),' Hz']);
    if (attitudeFlag == 2)
        disp(['MAGNETOMETER measurements at ',num2str(f_s), ' Hz']);
    else
        disp(['COMPASS measurements at ',num2str(f_s), ' Hz']);
    end
    disp('-------------------------------------------------------------------');
    disp('Simulating...');

end

end