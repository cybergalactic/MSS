% This function computes the 6-DOF generalized 1st-order wave forces, 
% tau_wave1, on a marine craft  using different wave spectra 
% (Modified Pierson-Moskowitz, JONSWAP, and Torsethaugen) and Response 
% Amplitude Operators (RAOs); see Fossen (2021, Chapters 10.2.1 and 10.2.4). 
% The real and imaginary parts of the RAO tables are interpolated in 
% frequency and varying wave directions to compute the RAO amplitudes and 
% phases. This approach avoids unwrapping problems and interpolation 
% issues in RAO phase angles.
displayMfileHeader('exWaveForceRAO.m');  % Print the header text

% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-15
% Revisions: 
%    2025-10-21 RAO interpolations and look-up tables are evaluated at a 
%               maximum rate of 10 Hz. 
%    2025-11-12 Added spheriod-shaped AUV with analytical RAO computations.
%    2026-05-14 Added spectrum parameters to GUI.

clearvars; close all;
clear waveForceRAO; % Clear persistent RAO tables
rng(1); % Set random generator seed to 1 when generating stochastic waves

[matFile, spectrumType, spreadingFlag, Hs, w0, beta_wave] = simOptions(); % GUI inputs
vessel = loadOrComputeVessel(matFile); % Load or compute vessel.forceRAO

% ------------------------------------------------------------------------------
% Simulation parameters
% ------------------------------------------------------------------------------
h = 0.02;                       % Time step [s]
T_final = 100;                  % Duration of the simulation [s]
T_initialTransient = 20;        % Remove initial transient [s]
RAO_update_period = 0.1;        % Compute RAO at 10 Hz

% COMMENT: Adding a spreading function involves summing waves from different 
% directions. Initially, these waves can interfere constructively, causing
% higher amplitudes. Hence, it is recommended to remove the initial
% respons by specifying: T_initialTransient >= 20 s.

% numFreqIntervals - Number of frequency intervals in wave spetrcum S(Omega)  
% numDirctions     - Number of wave directions in directional spectrum M(mu)
maxFreq = 3.0;                  % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 100;         % Number of wave frequency intervals (>50)
numDirections = 24;             % Number of wave directions (>15)

% ------------------------------------------------------------------------------
% Compute wave directional spectrum
% ------------------------------------------------------------------------------
spectrumParameters = [Hs, w0];
if strcmp(spectrumType ,'JONSWAP')
   gamma = 3.3; 
   spectrumParameters = [Hs, w0, gamma];
end

% Reshape vessel data to use 0 to maxFreq
if vessel.forceRAO.w(end) > maxFreq
    w_index = find(vessel.forceRAO.w > maxFreq, 1) - 1;
    vessel.forceRAO.w = vessel.forceRAO.w(1:w_index); % Frequency vector
    for DOF = 1:length(vessel.forceRAO.amp)
        vessel.forceRAO.amp{DOF} = vessel.forceRAO.amp{DOF}(1:w_index, :, :);
        vessel.forceRAO.phase{DOF} = vessel.forceRAO.phase{DOF}(1:w_index, :, :);
    end
end

% Calculate the wave spectrum power intensity S(Omega) for each frequency
omegaMax = vessel.forceRAO.w(end);  % Max frequency in RAO dataset
[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParameters, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

% ------------------------------------------------------------------------------
% MAIN LOOP
% ------------------------------------------------------------------------------
t = 0:h:T_final+T_initialTransient-1;  % Time vector
nextRAOtime = 0;                       % Next RAO update time
simdata = zeros(length(t),7);          % Pre-allocate table
for i = 1:length(t)

    U = 5;                             % Time-varying ship speed (m/s)
    psi = deg2rad(sin(0.1 * t(i)));    % Time-varying heading angle (rad)

    % 6-DOF generalized wave forces (compute RAO only every 0.1 second)
    if t(i) >= nextRAOtime
        [tau_wave1, waveElevation] = waveForceRAO(t(i), ...
            S_M, Amp, Omega, mu, vessel, U, psi, beta_wave, numFreqIntervals);

        nextRAOtime = nextRAOtime + RAO_update_period;
    end

    simdata(i,:) = [tau_wave1' waveElevation];

end

% ------------------------------------------------------------------------------
% PLOTS
% ------------------------------------------------------------------------------
figure(1); clf;

% Time-series
startIndex = max(1, floor(T_initialTransient / h) + 1);
t = t(startIndex:end) - t(startIndex);
tau_wave1 = simdata(startIndex:end, 1:6);
waveElevation = simdata(startIndex:end, 7);

% Plot the wave spectrum
subplot(211);
hold on;
if spreadingFlag
    % Plot the wave spectrum for the specific directions
    hold on;
    plot(Omega, S_M(:, floor(length(mu)/2)),'b->','Markersize',5,'LineWidth',1.5);
    plot(Omega, S_M(:, floor(length(mu)/4)),'k-o','Markersize',5, 'LineWidth',1.5);    
    plot(Omega, S_M(:, length(mu)),'g-', 'LineWidth',2.0);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))],'r-.', 'LineWidth', 1.5)
    legend('\mu = 0 deg', '\mu = 45 deg', '\mu = 90 deg',...
        ['\omega_0 = ', num2str(w0), ' rad/s']);
    hold off;
else
    hold on
    plot(Omega, S_M(:, 1), 'b-', 'LineWidth', 1.5);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))], 'r-.','LineWidth', 1.5)
    legend('S(\omega)', ['\omega_0 = ', num2str(w0), ' rad/s']);
    hold off
end
xlabel('Omega (rad/s)');
ylabel('m^2 s');
title([spectrumType, ' spectrum']);
grid on;

% Plot the wave elevation
subplot(212);
plot(t, waveElevation, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('m');
title(['Wave elevation for wave direction \beta_{wave} = ', ...
    num2str(rad2deg(beta_wave)), '°, H_s = ', num2str(Hs), ' m and ' ...
    '\omega_0 = ', num2str(w0), ' rad/s']);
grid on;

figure(2); clf;
% Plot the 6-DOF 1st-order wave forces
DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)',...
    'Roll (Nm)', 'Pitch (Nm)', 'Yaw (Nm)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, T_scale(DOF) * tau_wave1(:, DOF), 'LineWidth', 1.5);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF});
end

if ~isoctave
    sgtitle(['Generalized 1st-order wave-induced forces for \beta_{wave} = ' ...
        num2str(rad2deg(beta_wave)), '°, H_s = ', num2str(Hs), ' m and ' ...
        '\omega_0 = ', num2str(w0), ' rad/s'],'FontSize', 12);
end

% ------------------------------------------------------------------------------
% FUNCTIONS
% ------------------------------------------------------------------------------
function vessel = loadOrComputeVessel(matFile)
if strcmp(matFile, 'AUV_FUNCTION')
    % Example AUV parameters 
    a = 1.2;                   % Semi-major axis [m]
    b = 0.3;                   % Semi-minor axis [m]
    zn = 5.0;                  % Nominal submergence (m)

    vessel.main.g = 9.81;
    vessel = spheroidRAO(vessel,a,b,zn,0);
    disp('Generated vessel.forceRAO structure from function spheroidRAO().');

else
    load(which(matFile), 'vessel');
    disp(['Loaded the force RAO structure "vessel.forceRAO" (', matFile, ')'])
end
end

% ------------------------------------------------------------------------------

function [matFile, spectrumType, spreadingFlag, Hs, w0, beta_wave] = simOptions()

f = figure('Position', [400, 400, 450, 620], ...
    'Name', 'Simulation Options', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal');

% Add button group for selecting the ship
bg1 = uibuttongroup('Parent', f, ...
    'Position', [0.05 0.58 0.9 0.34], ...
    'Title', 'Select Marine Craft', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');
ship1 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Supply vessel', ...
    'Position', [10 90 150 30], ...
    'Tag', 'supply.mat');
ship2 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'S175 container ship', ...
    'Position', [10 60 150 30], ...
    'Tag', 's175.mat');
ship3 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Tanker', ...
    'Position', [10 30 150 30], ...
    'Tag', 'tanker.mat');
ship4 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'FPSO', ...
    'Position', [200 90 150 30], ...
    'Tag', 'FPSO.mat');
ship5 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Semisubmersible', ...
    'Position', [200 60 150 30], ...
    'Tag', 'semisub.mat');
ship6 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Spheriod-shaped AUV', ...
    'Position', [200 30 200 30], ...
    'Tag', 'AUV_FUNCTION');   

% Add button group for selecting the wave spectrum
bg2 = uibuttongroup('Parent', f, ...
    'Position', [0.05 0.33 0.9 0.23], ...
    'Title', 'Select Wave Spectrum', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');
spec1 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Modified Pierson-Moskowitz (PM)', ...
    'Position', [10 80 250 30], ...
    'Tag', 'Modified PM');
spec2 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'JONSWAP', ...
    'Position', [10 50 250 30], ...
    'Tag', 'JONSWAP', ...
    'Value', 1);
spec3 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Torsethaugen', ...
    'Position', [10 20 250 30], ...
    'Tag', 'Torsethaugen');

% Add checkbox for wave spreading
spreadCheckbox = uicontrol('Parent', f, ...
    'Style', 'checkbox', ...
    'FontSize', 11, ...
    'String', 'Enable Spreading Function', ...
    'Position', [35 180 250 30], ...
    'Value', 0);

% Add sea-state parameter inputs
uicontrol('Parent', f, ...
    'Style', 'text', ...
    'FontSize', 11, ...
    'String', 'Significant wave height Hs [m]', ...
    'HorizontalAlignment', 'left', ...
    'Position', [40 130 220 25]);

editHs = uicontrol('Parent', f, ...
    'Style', 'edit', ...
    'FontSize', 11, ...
    'String', '3', ...
    'Position', [300 128 80 28]);

uicontrol('Parent', f, ...
    'Style', 'text', ...
    'FontSize', 11, ...
    'String', 'Peak frequency w0 [rad/s]', ...
    'HorizontalAlignment', 'left', ...
    'Position', [40 95 220 25]);

editW0 = uicontrol('Parent', f, ...
    'Style', 'edit', ...
    'FontSize', 11, ...
    'String', '0.8', ...
    'Position', [300 93 80 28]);

uicontrol('Parent', f, ...
    'Style', 'text', ...
    'FontSize', 11, ...
    'String', 'Wave direction beta_wave [deg]', ...
    'HorizontalAlignment', 'left', ...
    'Position', [40 60 220 25]);

editBetaWave = uicontrol('Parent', f, ...
    'Style', 'edit', ...
    'FontSize', 11, ...
    'String', '140', ...
    'Position', [300 58 80 28]);

% Add OK button to confirm selections
uicontrol('Style', 'pushbutton', ...
    'String', 'OK', ...
    'FontSize', 12, ...
    'Position', [175 15 100 30], ...
    'Callback', @(src, evt) uiresume(f));


uiwait(f); % Wait for uiresume to be called on figure handle

Hs = str2double(get(editHs, 'String'));
w0 = str2double(get(editW0, 'String'));
beta_wave = deg2rad(str2double(get(editBetaWave, 'String')));

% Determine which ship was selected
selectedShip = findobj(bg1, 'Style', 'radiobutton', 'Value', 1);
matFile = get(selectedShip, 'Tag');

% Determine which wave spectrum was selected
selectedSpectrum = findobj(bg2, 'Style', 'radiobutton', 'Value', 1);
spectrumType = get(selectedSpectrum, 'Tag');

% Determine the state of the spreading checkbox
spreadingFlag = get(spreadCheckbox, 'Value');

if isnan(Hs) || Hs <= 0
    Hs = 3;
end
if isnan(w0) || w0 <= 0
    w0 = 0.8;
end
if isnan(beta_wave)
    beta_wave = deg2rad(140);
end

close(f); % Close the figure after obtaining the selections

end