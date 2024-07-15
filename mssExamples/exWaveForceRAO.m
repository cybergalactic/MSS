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

clear waveForceRAO; % Clear persistent RAO tables
clearvars;

[matFile, spectrumType, spreadingFlag] = simOptions(); % User inputs
load(which(matFile), 'vessel'); % Load vessel.forceRAO data structure
disp(['Loaded the force RAO structure "vessel.forceRAO" (', matFile, ')'])

% Simulation parameters
h = 0.1;                        % Time step (s)
T_final = 200;                  % Duration of the simulation (s)
T_initialTransient = 20;        % Remove initial transient (s)

% COMMENT: Adding a spreading function involves summing waves from different 
% directions. Initially, these waves can interfere constructively, causing
% higher amplitudes. Hence, it is recommended to remove the initial
% respons by specifying: T_initialTransient >= 20 s.

% numFreqIntervals - Number of frequency intervals in wave spetrcum S(Omega)  
% numDirctions     - Number of wave directions in directional spectrum M(mu)
numFreqIntervals = 100;         % Number of wave frequency intervals (>100)
numDirections = 24;             % Number of wave directions (>15)

% Sea state 
Hs = 10;                        % Significant wave height (m)
Tz = 10;                        % Zero-crossing period (s)
    
% Wave direction relative bow, 0 deg for following sea, 180 deg for head sea
beta_wave = deg2rad(140);

% Calculate the wave spectrum power intensity S(Omega) for each frequency
T0 = Tz / 0.710; % Wave spectrum modal (peak) period (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;  % Wave spectrum modal (peak) frequency
spectrumParameters = [Hs, w0];
if strcmp(spectrumType ,'JONSWAP')
   gamma = 3.3; 
   spectrumParameters = [Hs, w0, gamma];
end

omegaMax = vessel.forceRAO.w(end);  % Max frequency in RAO dataset

[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParameters, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% MAIN LOOP
t = 0:h:T_final+T_initialTransient-1;  % Time vector
simdata = zeros(length(t),7);          % Pre-allocate table
for i = 1:length(t)

    U = 5;                             % Time-varying ship speed (m/s)
    psi = deg2rad(sin(0.1 * t(i)));    % Time-varying heading angle (rad)

    % 6-DOF generalized wave forces
    [tau_wave1, waveElevation] = waveForceRAO(t(i), ...
        S_M, Amp, Omega, mu, vessel, U, psi, beta_wave, numFreqIntervals);

    simdata(i,:) = [tau_wave1' waveElevation];

end

%% PLOTS
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
    plot(Omega, S_M(:, floor(length(mu)/2)), 'LineWidth', 2);
    plot(Omega, S_M(:, floor(length(mu)/4)), 'LineWidth', 2);    
    plot(Omega, S_M(:, length(mu)), 'LineWidth', 2);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))], 'LineWidth', 2)
    legend('\mu = 0 deg', '\mu = 45 deg', '\mu = 90 deg',...
        ['w_0 = ', num2str(w0), ' rad/s']);
    hold off;
else
    hold on
    plot(Omega, S_M(:, 1), 'b-', 'LineWidth', 2);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))], 'LineWidth', 2)
    legend('S(Omega)', ['w_0 = ', num2str(w0), ' rad/s']);
    hold off
end
xlabel('Omega (rad/s)');
ylabel('m^2 s');
title([spectrumType, ' spectrum']);
grid on;

% Plot the wave elevation
subplot(212);
plot(t, waveElevation, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('m');
title(['Wave Elevation for wave direction \beta_{wave} = ', ...
    num2str(rad2deg(beta_wave)), '°']);
grid on;

figure(2); clf;
% Plot the 6-DOF 1st-order wave forces
DOF_txt = {'Surge (m)', 'Sway (m)', 'Heave (m)',...
    'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, T_scale(DOF) * tau_wave1(:, DOF), 'LineWidth', 2);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF});
end

if ~isoctave
    sgtitle(['Generalized 1st-order Wave Forces for \beta_{wave} = ' ...
        num2str(rad2deg(beta_wave)), '°']);
end


%% FUNCTIONS
function [matFile, spectrumType, spreadingFlag] = simOptions()

f = figure('Position', [400, 400, 400, 450], ...
    'Name', 'Simulation Options', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal');

% Add button group for selecting the ship
bg1 = uibuttongroup('Parent', f, ...
    'Position', [0.05 0.55 0.9 0.35], ...
    'Title', 'Select Ship', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');
ship1 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Supply vessel', ...
    'Position', [10 80 150 30], ...
    'Tag', 'supply.mat');
ship2 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'S175 container ship', ...
    'Position', [10 50 150 30], ...
    'Tag', 's175.mat');
ship3 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Tanker', ...
    'Position', [10 20 150 30], ...
    'Tag', 'tanker.mat');
ship4 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'FPSO', ...
    'Position', [200 80 150 30], ...
    'Tag', 'FPSO.mat');
ship5 = uicontrol(bg1, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Semisubmersible', ...
    'Position', [200 50 150 30], ...
    'Tag', 'semisub.mat');

% Add button group for selecting the wave spectrum
bg2 = uibuttongroup('Parent', f, ...
    'Position', [0.05 0.3 0.9 0.25], ...
    'Title', 'Select Wave Spectrum', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');
spec1 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Modified Pierson-Moskowitz (PM)', ...
    'Position', [10 60 250 30], ...
    'Tag', 'Modified PM');
spec2 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'JONSWAP', ...
    'Position', [10 30 250 30], ...
    'Tag', 'JONSWAP');
spec3 = uicontrol(bg2, ...
    'Style', 'radiobutton', ...
    'FontSize', 11, ...
    'String', 'Torsethaugen', ...
    'Position', [10 0 250 30], ...
    'Tag', 'Torsethaugen');

% Add checkbox for wave spreading
spreadCheckbox = uicontrol('Parent', f, ...
    'Style', 'checkbox', ...
    'FontSize', 12, ...
    'String', 'Enable Spreading Function', ...
    'Position', [100 100 200 30], ...
    'Value', 0); % Default is unchecked

% Add OK button to confirm selections
uicontrol('Style', 'pushbutton', ...
    'String', 'OK', ...
    'FontSize', 12, ...
    'Position', [150 50 100 40], ...
    'Callback', @(src, evt) uiresume(f));

uiwait(f); % Wait for uiresume to be called on figure handle

% Determine which ship was selected
selectedShip = findobj(bg1, 'Style', 'radiobutton', 'Value', 1);
matFile = get(selectedShip, 'Tag');

% Determine which wave spectrum was selected
selectedSpectrum = findobj(bg2, 'Style', 'radiobutton', 'Value', 1);
spectrumType = get(selectedSpectrum, 'Tag');

% Determine the state of the spreading checkbox
spreadingFlag = get(spreadCheckbox, 'Value');

close(f); % Close the figure after obtaining the selections

end