function SIMaidedINSeuler()
% SIMaidedINSeuler is compatible with MATLAB and GNU Octave (www.octave.org).
%
% This function simulates two error-state Kalman filter (ESKF) architectures 
% for an inertial navigation system (INS):
%
%   1. A 15-state ESKF in which attitude and ARS bias are estimated using
%      IMU, compass, position, and optionally velocity measurements.
%
%   2. A 9-state ESKF in which attitude is supplied by an external
%      attitude and heading reference system (AHRS). The ESKF estimates
%      position, velocity, and accelerometer bias using position and
%      optionally velocity measurements.
%
% Attitude is parameterized using Euler angles (Fossen, 2027, Chapter 14).
%
% The aiding frequency f_slow can be selected independently of the high-rate 
% IMU frequency f_fast, subject to f_slow <= f_fast.
%
% Dependencies:
%   ins_euler.m
%       Feedback ESKF for an INS aided by compass and position  measurements. 
%       Velocity aiding is optional.
%
%   ins_ahrs.m
%       Feedback ESKF for an INS aided by an external AHRS and position
%       measurements. Velocity aiding is optional.
%
%   insSignal.m
%       INS signal generator.
%
%   magneticField.m
%       Magnetic-field reference vectors and latitude for selected cities.
%
% References:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and
%   Motion Control, 3rd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
% Author: Thor I. Fossen
% Date: 2021-04-26
% Revisions:
%   2024-08-20: Using the updated insSignal.m generator.
%   2024-11-02: Improved logic for slow position data.
%   2026-07-14: Added support for both the 15-state compass-aided ESKF
%               and the 9-state AHRS-assisted ESKF.

% ==============================================================================
% Simulation parameters
% ==============================================================================
T_final = 100;        % Final simulation time (s)
f_fast  = 1000;       % IMU and filter propagation frequency (Hz)
f_slow  = 5;          % Position/velocity aiding frequency (Hz)

h      = 1 / f_fast;  % High-rate IMU and ESKF sampling time
h_slow = 1 / f_slow;  % Low-rate aiding sampling time

testSignalNo = 1;     % INS test signal - 1: constant bias, 2: time-varying bias

% ==============================================================================
% Initialization of the ESKF 
% ==============================================================================
p_0 = 1.0; % Initial covariance matrix: P_prd = p_0 * I_nxn

% Measurement standard deviations
sigma_pos = 0.05;       % Position [m]
sigma_vel = 0.01;       % Velocity [m/s]
sigma_psi = deg2rad(1); % Compass heading [rad]
sigma_g   = 0.1;        % Normalized gravity-vector residual

% Process noise
q_f     = 1e-3;         % Specific-force process noise
q_b_acc = 1e-5;         % Accelerometer bias process noise
q_w     = 1e-3;         % Angular-rate process noise
q_b_ars = 1e-5;         % ARS bias process noise

% Bias time constants
T_acc = 300;            % Acceleration bias time constant [s]
T_ars = 300;            % Angular rate bias time constant [s]

% ==============================================================================
% Initialization of the INS signal generator
% ==============================================================================
[m_ref, ~, mu, cityName] = magneticField(1);

% True IMU biases used by the signal generator
b_acc = [0.1  0.3  -0.1]';
b_ars = [0.05 0.1  -0.05]';

% Signal-generator state:
% [position; velocity; accelerometer bias; Euler angles; ARS bias]
x_true = [zeros(1,6), b_acc', zeros(1,3), b_ars']';

% Select attitude source and optional velocity aiding
[attitudeMethod, aidingMethod] = displayMethod(cityName, f_fast, f_slow);

% ==============================================================================
% Initialization of the ESKF covariance matrices
% ==============================================================================
switch attitudeMethod

    case 'compass'
        % Initialize the 15-state ESKF covariance matrices
        % delta_x = [delta_p; delta_v; delta_b_acc; delta_theta; delta_b_ars]
        P_prd_compass = p_0 * eye(15);

        Qd = diag([ ...
            q_f q_f q_f ...
            q_b_acc q_b_acc q_b_acc ...
            q_w q_w q_w ...
            q_b_ars q_b_ars q_b_ars ]);

        % Measurement noise:
        switch aidingMethod

            case 'position'
                % [position; normalized specific force; compass heading]
                Rd = diag([ ...
                    sigma_pos^2 sigma_pos^2 sigma_pos^2 ...
                    sigma_g^2   sigma_g^2   sigma_g^2 ...
                    sigma_psi^2 ]);

            case 'position_velocity'
                % [position; velocity; normalized specific force; compass heading]
                Rd = diag([ ...
                    sigma_pos^2 sigma_pos^2 sigma_pos^2 ...
                    sigma_vel^2 sigma_vel^2 sigma_vel^2 ...
                    sigma_g^2   sigma_g^2   sigma_g^2 ...
                    sigma_psi^2 ]);

            otherwise
                error('Unknown aiding method: %s', aidingMethod)
        end

    case 'ahrs'

        % Initialize the 9-state ESKF covariance matrices
        % delta_x = [delta_p; delta_v; delta_b_acc]
        P_prd_ahrs = p_0 * eye(9);

        % Process noise
        Qd = diag([ ...
            q_f q_f q_f ...
            q_b_acc q_b_acc q_b_acc ]);

        switch aidingMethod
            case 'position'

                % Measurement noise: position
                Rd = sigma_pos^2 * eye(3);

            case 'position_velocity'

                % Measurement noise: position and velocity
                Rd = diag([ ...
                    sigma_pos^2 sigma_pos^2 sigma_pos^2 ...
                    sigma_vel^2 sigma_vel^2 sigma_vel^2 ]);

            otherwise
                error('Unknown aiding method: %s', aidingMethod)

        end

    otherwise
        error('Unknown attitude method: %s', attitudeMethod)
end

% ==============================================================================
% Initialization of the INS states
% ==============================================================================
p_ins     = zeros(3,1);
v_ins     = zeros(3,1);
b_acc_ins = zeros(3,1);
theta_ins = zeros(3,1);
b_ars_ins = zeros(3,1);

switch attitudeMethod

    case 'compass'
        % State propagated by ins_euler.m (15 states)
        x_ins_compass = [ ...
            p_ins;
            v_ins;
            b_acc_ins;
            theta_ins;
            b_ars_ins];

    case 'ahrs'
        % State propagated by ins_ahrs.m (9 states)
        x_ins_ahrs = [ ...
            p_ins;
            v_ins;
            b_acc_ins];
end

% ==============================================================================
% Multirate scheduling
% ==============================================================================
% At most one slow measurement is processed per fast time step.
if f_slow > f_fast
    error('The aiding frequency f_slow must satisfy f_slow <= f_fast.');
end

slowIndex = 0; % Zero-based index of the next nominal slow measurement

% Tolerance for floating-point comparisons of coincident sample times
tol = 10 * eps(max(1, T_final));

% ==============================================================================
% Time and data initialization
% ==============================================================================
t = 0:h:T_final;
nTimeSteps = length(t);
maxSlowSamples = floor(T_final / h_slow) + 1;

trueData = zeros(nTimeSteps, 15); % True states: % [p; v; b_acc; theta; b_ars]
navEstimateData = zeros(nTimeSteps, 9); % Navigation estimates: % [p; v; b_acc]

% Attitude data. In compass mode these are ESKF estimates; in AHRS mode
% they are the external AHRS measurements used by the navigation filter
attitudeData = zeros(nTimeSteps, 3);

% ARS-bias estimates exist only for the compass-aided 15-state ESKF
arsBiasEstimateData = nan(nTimeSteps, 3);

% Slow position measurements
positionData = zeros(maxSlowSamples, 4);
positionIndex = 0;

% ==============================================================================
%% MAIN LOOP
% ==============================================================================
for i = 1:nTimeSteps

    % INS signal generator
    [x_true, f_imu, w_imu] = insSignal(x_true, h, t(i), mu, m_ref, testSignalNo);

    % Compass and AHRS outputs
    y_psi  = x_true(12);
    y_ahrs = x_true(10:12);

    % Determine whether a new slow aiding measurement is available
    newSlowMeasurement = t(i) + tol >= slowIndex * h_slow;

    if newSlowMeasurement

        slowIndex = slowIndex + 1;
        positionIndex = positionIndex + 1;

        y_pos = x_true(1:3) + 0.05 * randn(3,1);
        y_vel = x_true(4:6) + 0.01 * randn(3,1);

        positionData(positionIndex,:) = [t(i), y_pos'];

        switch attitudeMethod
            case 'compass'

                switch aidingMethod

                    case 'position'
                        [x_ins_compass, P_prd_compass] = ins_euler( ...
                            x_ins_compass,P_prd_compass,mu,h,Qd,Rd, ...
                            T_acc,T_ars,[f_imu' w_imu'],y_psi,y_pos);

                    case 'position_velocity'

                        [x_ins_compass, P_prd_compass] = ins_euler( ...
                            x_ins_compass,P_prd_compass,mu,h,Qd,Rd, ...
                            T_acc,T_ars,[f_imu' w_imu'],y_psi,y_pos,y_vel);
                end

            case 'ahrs'

                switch aidingMethod
                    case 'position'

                        [x_ins_ahrs, P_prd_ahrs] = ins_ahrs( ...
                            x_ins_ahrs,P_prd_ahrs,mu,h,Qd,Rd, ...
                            T_acc,f_imu,y_ahrs,y_pos);

                    case 'position_velocity'

                        [x_ins_ahrs, P_prd_ahrs] = ins_ahrs( ...
                            x_ins_ahrs,P_prd_ahrs,mu,h,Qd,Rd, ...
                            T_acc,f_imu,y_ahrs,y_pos,y_vel);
                end

        end

    else

        % No new low-rate position or velocity aiding measurement
        switch attitudeMethod
            case 'compass'

                [x_ins_compass, P_prd_compass] = ins_euler( ...
                    x_ins_compass,P_prd_compass,mu,h,Qd,Rd, ...
                    T_acc,T_ars,[f_imu' w_imu'],y_psi);

            case 'ahrs'

                [x_ins_ahrs, P_prd_ahrs] = ins_ahrs( ...
                    x_ins_ahrs,P_prd_ahrs,mu,h,Qd,Rd,T_acc,f_imu,y_ahrs);
        end

    end

    % Store simulation data
    trueData(i,:) = x_true';

    switch attitudeMethod

        case 'compass'
            navEstimateData(i,:)      = x_ins_compass(1:9)';
            attitudeData(i,:)         = x_ins_compass(10:12)';
            arsBiasEstimateData(i,:)  = x_ins_compass(13:15)';
        
        case 'ahrs'
            navEstimateData(i,:) = x_ins_ahrs';
            attitudeData(i,:) = y_ahrs';
    end

end

% Remove unused preallocated rows
positionData = positionData(1:positionIndex,:);

% ==============================================================================
%%  PLOTS
% ==============================================================================
scrSz = get(0, 'ScreenSize');

legendSize = 10;
colors = {'b', 'g', 'k'};

xTrue = trueData;
xNav  = navEstimateData;

t_m = positionData(:,1);
y_m = positionData(:,2:4);

% ==============================================================================
% Figure 1: Translational navigation states
% ==============================================================================
figure(1);
clf

if ~isoctave
    set(gcf, 'Position', [1, 1, 0.4 * scrSz(3), scrSz(4)]);
end

% ------------------------------------------------------------------------------
% Position
% ------------------------------------------------------------------------------
subplot(3,1,1)

hMeas = plot(t_m, y_m, 'xr');
hold on

hX = plot(t, xNav(:,1), colors{1});
hY = plot(t, xNav(:,2), colors{2});
hZ = plot(t, xNav(:,3), colors{3});

hAll = [hX, hY, hZ, hMeas(1)];

hold off
xlabel('Time [s]')
title('Position [m]')
grid on

labels = { ...
    ['Estimate x_N at ', num2str(f_fast), ' Hz'], ...
    ['Estimate y_E at ', num2str(f_fast), ' Hz'], ...
    ['Estimate z_D at ', num2str(f_fast), ' Hz'], ...
    ['Position measurements at ', num2str(f_slow), ' Hz']};

legend(hAll, labels)

% ------------------------------------------------------------------------------
% Velocity
% ------------------------------------------------------------------------------
subplot(3,1,2)

hTrue = plot(t, xTrue(:,4:6), 'r');
hold on

hX = plot(t, xNav(:,4), colors{1});
hY = plot(t, xNav(:,5), colors{2});
hZ = plot(t, xNav(:,6), colors{3});

hAll = [hX, hY, hZ, hTrue(1)];

hold off
xlabel('Time [s]')
title('Velocity [m/s]')
grid on

labels = { ...
    ['Estimate v_N at ', num2str(f_fast), ' Hz'], ...
    ['Estimate v_E at ', num2str(f_fast), ' Hz'], ...
    ['Estimate v_D at ', num2str(f_fast), ' Hz'], ...
    ['True velocity at ', num2str(f_fast), ' Hz']};

legend(hAll, labels)

% ------------------------------------------------------------------------------
% Accelerometer bias
% ------------------------------------------------------------------------------
subplot(3,1,3)

hTrue = plot(t, xTrue(:,7:9), 'r');
hold on

hX = plot(t, xNav(:,7), colors{1});
hY = plot(t, xNav(:,8), colors{2});
hZ = plot(t, xNav(:,9), colors{3});

hAll = [hX, hY, hZ, hTrue(1)];

hold off
xlabel('Time [s]')
title('Accelerometer bias [m/s^2]')
grid on

labels = { ...
    ['Estimate b_{x,acc} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate b_{y,acc} at ', num2str(f_fast), ' Hz'], ...
    ['Estimate b_{z,acc} at ', num2str(f_fast), ' Hz'], ...
    ['True accelerometer bias at ', num2str(f_fast), ' Hz']};

legend(hAll, labels)

set(findall(gcf, 'Type', 'line'),   'LineWidth', 1.5)
set(findall(gcf, 'Type', 'text'),   'FontSize', 12)
set(findall(gcf, 'Type', 'legend'), 'FontSize', legendSize)

% ==============================================================================
% Figure 2: Attitude and ARS bias
% ==============================================================================
figure(2);
clf

if ~isoctave
    set(gcf, 'Position', ...
        [0.4 * scrSz(3), 1, 0.4 * scrSz(3), scrSz(4)]);
end

switch attitudeMethod
    case 'compass'

        % --------------------------------------------------------------------------
        % Euler-angle estimates from the 15-state ESKF
        % --------------------------------------------------------------------------

        subplot(2,1,1)

        hTrue = plot(t, rad2deg(xTrue(:,10:12)), 'r');
        hold on

        hPhi   = plot(t, rad2deg(attitudeData(:,1)), colors{1});
        hTheta = plot(t, rad2deg(attitudeData(:,2)), colors{2});
        hPsi   = plot(t, rad2deg(attitudeData(:,3)), colors{3});

        hAll = [hPhi, hTheta, hPsi, hTrue(1)];

        hold off
        xlabel('Time [s]')
        title('Euler angles [deg]')
        grid on

        labels = { ...
            ['Estimate \phi at ', num2str(f_fast), ' Hz'], ...
            ['Estimate \theta at ', num2str(f_fast), ' Hz'], ...
            ['Estimate \psi at ', num2str(f_fast), ' Hz'], ...
            ['True Euler angles at ', num2str(f_fast), ' Hz']};

        legend(hAll, labels)

        % --------------------------------------------------------------------------
        % ARS-bias estimates from the 15-state ESKF
        % --------------------------------------------------------------------------
        subplot(2,1,2)

        hTrue = plot(t, rad2deg(xTrue(:,13:15)), 'r');
        hold on

        hX = plot(t, rad2deg(arsBiasEstimateData(:,1)), colors{1});
        hY = plot(t, rad2deg(arsBiasEstimateData(:,2)), colors{2});
        hZ = plot(t, rad2deg(arsBiasEstimateData(:,3)), colors{3});

        hAll = [hX, hY, hZ, hTrue(1)];

        hold off
        xlabel('Time [s]')
        title('Angular-rate bias [deg/s]')
        grid on

        labels = { ...
            ['Estimate b_{x,ars} at ', num2str(f_fast), ' Hz'], ...
            ['Estimate b_{y,ars} at ', num2str(f_fast), ' Hz'], ...
            ['Estimate b_{z,ars} at ', num2str(f_fast), ' Hz'], ...
            ['True ARS bias at ', num2str(f_fast), ' Hz']};

        legend(hAll, labels)

    case 'ahrs'

        % --------------------------------------------------------------------------
        % Attitude supplied by the external AHRS
        % --------------------------------------------------------------------------
        hTrue = plot(t, rad2deg(xTrue(:,10:12)), 'r');
        hold on

        hPhi   = plot(t, rad2deg(attitudeData(:,1)), colors{1});
        hTheta = plot(t, rad2deg(attitudeData(:,2)), colors{2});
        hPsi   = plot(t, rad2deg(attitudeData(:,3)), colors{3});

        hAll = [hPhi, hTheta, hPsi, hTrue(1)];

        hold off
        xlabel('Time [s]')
        title('AHRS Euler angles [deg]')
        grid on

        labels = { ...
            ['AHRS \phi at ', num2str(f_fast), ' Hz'], ...
            ['AHRS \theta at ', num2str(f_fast), ' Hz'], ...
            ['AHRS \psi at ', num2str(f_fast), ' Hz'], ...
            ['True Euler angles at ', num2str(f_fast), ' Hz']};

        legend(hAll, labels)

end

set(findall(gcf, 'Type', 'line'),   'LineWidth', 1.5)
set(findall(gcf, 'Type', 'text'),   'FontSize', 12)
set(findall(gcf, 'Type', 'legend'), 'FontSize', legendSize)


% ==============================================================================
% Radio buttons, flags, and display
% ==============================================================================
function [attitudeMethod, aidingMethod] = displayMethod( ...
    cityName, f_fast, f_slow)

    f = figure( ...
        'Position',    [400, 400, 500, 320], ...
        'Name',        'Strapdown Aided INS', ...
        'MenuBar',     'none', ...
        'NumberTitle', 'off', ...
        'WindowStyle', 'modal');

    % --------------------------------------------------------------------------
    % Attitude source
    % --------------------------------------------------------------------------
    bg1 = uibuttongroup( ...
        'Parent',     f, ...
        'Position',   [0.02, 0.64, 0.96, 0.32], ...
        'Title',      'Attitude Source', ...
        'FontSize',   14, ...
        'FontWeight', 'bold');

    radioCompass = uicontrol( ...
        bg1, ...
        'Style',    'radiobutton', ...
        'FontSize', 13, ...
        'String',   'Compass-aided 15-state ESKF', ...
        'Position', [10, 45, 470, 30], ...
        'Tag',      'compass');

    radioAHRS = uicontrol( ...
        bg1, ...
        'Style',    'radiobutton', ...
        'FontSize', 13, ...
        'String',   'External AHRS with 9-state ESKF', ...
        'Position', [10, 12, 470, 30], ...
        'Tag',      'ahrs');

    set(radioCompass, 'Value', 1);

    % --------------------------------------------------------------------------
    % Aiding measurements
    % --------------------------------------------------------------------------
    bg2 = uibuttongroup( ...
        'Parent',     f, ...
        'Position',   [0.02, 0.32, 0.96, 0.28], ...
        'Title',      'Aiding Measurements', ...
        'FontSize',   14, ...
        'FontWeight', 'bold');

    radioPosition = uicontrol( ...
        bg2, ...
        'Style',    'radiobutton', ...
        'FontSize', 13, ...
        'String',   'Position aiding', ...
        'Position', [10, 35, 470, 30], ...
        'Tag',      'position');

    radioPositionVelocity = uicontrol( ...
        bg2, ...
        'Style',    'radiobutton', ...
        'FontSize', 13, ...
        'String',   'Position and velocity aiding', ...
        'Position', [10, 5, 470, 30], ...
        'Tag',      'position_velocity');

    set(radioPosition, 'Value', 1);

    % --------------------------------------------------------------------------
    % Confirmation button
    % --------------------------------------------------------------------------
    uicontrol( ...
        'Style',    'pushbutton', ...
        'String',   'OK', ...
        'FontSize', 13, ...
        'Position', [20, 25, 100, 40], ...
        'Callback', @(src, event) uiresume(f));

    uiwait(f);

    % Determine the selected attitude source
    if get(radioCompass, 'Value') == 1
        attitudeMethod = get(radioCompass, 'Tag');
    else
        attitudeMethod = get(radioAHRS, 'Tag');
    end

    % Determine the selected aiding measurements
    if get(radioPosition, 'Value') == 1
        aidingMethod = get(radioPosition, 'Tag');
    else
        aidingMethod = get(radioPositionVelocity, 'Tag');
    end

    close(f);

    % --------------------------------------------------------------------------
    % Display simulation configuration
    % --------------------------------------------------------------------------
    disp('-------------------------------------------------------------------')
    disp('MSS toolbox: Error-state feedback Kalman filter')
    disp('Attitude parameterization: Euler angles')

    switch aidingMethod
        case 'position'
            disp([ ...
                'INS aided by position measurements at ', ...
                num2str(f_slow), ' Hz'])

        case 'position_velocity'
            disp([ ...
                'INS aided by position and velocity measurements at ', ...
                num2str(f_slow), ' Hz'])
    end

    disp([ ...
        'IMU measurements and ESKF propagation at ', ...
        num2str(f_fast), ' Hz'])

    switch attitudeMethod
        case 'compass'
            disp('Architecture: 15-state ESKF with compass aiding')
            disp([ ...
                'Compass measurements at ', ...
                num2str(f_fast), ' Hz'])

        case 'ahrs'
            disp('Architecture: 9-state ESKF using an external AHRS')
            disp([ ...
                'Three-axis AHRS measurements at ', ...
                num2str(f_fast), ' Hz'])
    end

    disp([ ...
        'Magnetic-field reference vector for ', ...
        cityName, ' (>> type magneticField)'])

    disp('-------------------------------------------------------------------')
    disp('Simulating...')
end

end