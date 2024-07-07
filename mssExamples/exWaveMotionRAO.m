function exWaveMotionRAO()
% exWaveMotionRAO is compatibel with MATLAB and GNU Octave (www.octave.org).
% This function computes the wave elevation and the wave-frequency (WF) 
% motion, eta_w, on a ship using different wave spectra (Modified Pierson-
% Moskowitz, JONSWAP, Torsethaugen) and Response Amplitude Operators (RAOs) 
% (Fossen 2021, Chapters 10.2.1 and 10.2.3). The total motion is:
%
%   y = eta + eta_w
%
% where eta is the low-frequency (LF) ship model component. The real and 
% imaginary parts of the RAO tables are interpolated in frequency and 
% varying wave directions to compute the RAO amplitudes and phases. This 
% approach avoids unwrapping problems and interpolation issues in RAO 
% phase angles.
displayMfileHeader('exWaveMotionRAO.m');  % Print the header text

% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-06
% Revisions: 

[matFile, spectrumType, spreadingFlag] = simOptions(); % User inputs
load(which(matFile), 'vessel'); % Load vessel.motionRAO data structure
disp(['Loaded the motion RAO structure "vessel.motionRAO" (', matFile, ')'])

%% Simulation parameters
Tfinal = 200;                   % Duration of the simulation (s)
h = 0.1;                        % Time step (s)
numFreqIntervals = 100;         % Number of wave frequency intervals (>100)
g = 9.81;                       % Acceleration of gravity (m/s^2)

% Sea state - JONSWAP spectrum with DNV (2007) correction
Hs = 10;                        % Significant wave height (m)
Tz = 10;                        % Zero-crossing period (s)
    
% Wave direction relative bow, 0 deg for following sea, 180 deg for head sea
beta_wave = deg2rad(140);

% Fixed seed for reproducibility of random phase angles
seed = 12345; 
rng(seed);

freqs = vessel.motionRAO.w;             % RAO wave frequencies
raoAngles = vessel.headings;            % RAO wave direction angles
numAngles = length(raoAngles);          
omegaMax = freqs(end);

% Calculate the wave spectrum power intensity S(Omega) for each frequency
T0 = Tz / 0.710; % Wave spectrum modal (peak) period (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;  % Wave spectrum modal (peak) frequency
spectrumParameters = [Hs, w0];
if strcmp(spectrumType ,'JONSWAP')
   gamma = 3.3; 
   spectrumParameters = [Hs, w0, gamma];
end

[S_M, Omega, Amp, ~, ~, mu] = waveSpectrum(spectrumType, ...
    spectrumParameters, numFreqIntervals, omegaMax, spreadingFlag);

% Convert RAO amplitudes and phases to real and imaginary components
RAO_re = cell(6, numAngles);
RAO_im = cell(6, numAngles);
for DOF = 1:6
    for j = 1:numAngles
        % Extract amp and phase for zero speed, index 1
        RAO_phase = vessel.motionRAO.phase{DOF}(:,:,1);
        RAO_amp = vessel.motionRAO.amp{DOF}(:,:,1);

        % Calculate the real and imaginary parts
        RAO_re{DOF, j} = RAO_amp .* cos(RAO_phase);
        RAO_im{DOF, j} = RAO_amp .* sin(RAO_phase);
    end
end

% Interpolate Re and Im parts of RAO to be valid for all Omega values
% Repeat this for all wave directions k = 1:numAngles
F_re_values_all = cell(1, 6); 
F_im_values_all = cell(1, 6);
for DOF = 1:6
    F_re_values = zeros(numAngles, numFreqIntervals);
    F_im_values = zeros(numAngles, numFreqIntervals);

    for k = 1:numAngles
        F_re_values(k, :) = interp1(freqs(:), RAO_re{DOF, k}(:, k), ...
            Omega, 'linear', 'extrap');
        F_im_values(k, :) = interp1(freqs(:), RAO_im{DOF, k}(:, k), ...
            Omega, 'linear', 'extrap');
    end

    % Store the interpolated values for this DOF
    F_re_values_all{DOF} = F_re_values;
    F_im_values_all{DOF} = F_im_values;
end

%% MAIN LOOP
t = 0:h:Tfinal; % Time vector
eta_WF = zeros(length(t), 6); 
waveElevation = zeros(length(t), 1);
randomPhases = 2 * pi * rand(numFreqIntervals, 1);

for t_i = 1:length(t)

    U = 5;                             % Time-varying ship speed (m/s)
    psi = deg2rad(sin(0.1 * t(t_i)));  % Time-varying heading angle (rad)

    beta_relative = beta_wave - psi;  % Wave direction relative ship
    
    % Vector of spreading angles, scalar for M = 1 and one direction mu = 0
    angle_spreading = mod(beta_relative - mu, 2*pi); % Wrap to 0 to 2*pi 

    % Encounter frequency for all frequencies and directions
    Omega_e = abs(Omega - (Omega.^2 / g) * U .* cos(angle_spreading'));

    % Compute the wave elevation (Fossen 2021, Eq. 10.83)
    waveElevation(t_i) = 0; % Initialize waveElevation for this time step

    if size(S_M, 2) == 1 % No spreading function
        integrand = Amp .* cos(Omega_e .* t(t_i) + randomPhases);
    else % Directional spectrum , apply trapz for summation over directions mu
        integrand = zeros(numFreqIntervals, 1);
        for w_i = 1:numFreqIntervals
            integrand(w_i) = trapz(mu, Amp(w_i, :) .* ...
                 cos(Omega_e(w_i, :) .* t(t_i) + randomPhases(w_i)));
        end
    end
    % Apply trapz for summation over frequenies Omega
    waveElevation(t_i) = trapz(Omega, integrand);

    % Compute the complex RAOs as a function of frequency and wave direction
    RAO_complex = cell(1, 6); % Initialize cell arrays
    for DOF = 1:6

        % Retrieve stored interpolated values
        F_re_values = F_re_values_all{DOF};
        F_im_values = F_im_values_all{DOF};

        % Initialize interpolation results
        F_re_dir_interp = zeros(length(mu), numFreqIntervals);
        F_im_dir_interp = zeros(length(mu), numFreqIntervals);

        % Interpolate Re and Im parts of RAO for wave directions 0 to 2*pi
        for k = 1:length(mu)
            F_re_dir_interp(k, :) = interp1(raoAngles, F_re_values, ...
                angle_spreading(k), 'linear', 'extrap');
            F_im_dir_interp(k, :) = interp1(raoAngles, F_im_values, ...
                angle_spreading(k), 'linear', 'extrap');
        end

        % Combine real and imaginary parts to form the complex RAO
        RAO_complex{DOF} = F_re_dir_interp + 1i * F_im_dir_interp;

        % Compute the 6-DOF wave-frequency motions (Fossen 2021, Eq. 10.105).
        integrand = zeros(numFreqIntervals, 1);

        if size(S_M, 2) == 1 % No spreading function/directional spectrum

            for w_i = 1:numFreqIntervals
                integrand(w_i) = abs(RAO_complex{DOF}(:, w_i))' ...
                    .* Amp(w_i, :) .* cos(Omega_e(w_i, :) .* t(t_i) + ...
                    angle(RAO_complex{DOF}(:, w_i))' + randomPhases(w_i));
            end   

        else % Directional spectrum

            for w_i = 1:numFreqIntervals
                % Apply trapezoidal rule for summation over directions first
                integrand(w_i) = trapz(mu, ...
                    abs(RAO_complex{DOF}(:, w_i))' .* ...
                    Amp(w_i, :) .* cos(Omega_e(w_i, :) .* t(t_i) + ...
                    angle(RAO_complex{DOF}(:, w_i))' + randomPhases(w_i)));
            end

        end

        % Apply trapezoidal rule for summation over frequencies
        eta_WF(t_i, DOF) = trapz(Omega, integrand);

    end
end

%% PLOTS
figure(1); clf;
% Plot the wave spectrum
subplot(311);
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

% Plot the wave amplitude
subplot(312);
if spreadingFlag
    % Plot the wave spectrum for the specific directions
    hold on;
    plot(Omega, Amp(:, floor(length(mu)/2)), 'LineWidth', 2);
    plot(Omega, Amp(:, floor(length(mu)/4)), 'LineWidth', 2);    
    plot(Omega, Amp(:, length(mu)), 'LineWidth', 2);
    legend('\mu = 0 deg', '\mu = 45 deg', '\mu = 90 deg');
    hold off
else
    plot(Omega, Amp(:, 1), 'b-', 'LineWidth', 2);
end
xlabel('Omega (rad/s)');
ylabel('m');
title('Wave Amplitudes');
grid on;

% Plot the wave elevation
subplot(313);
plot(t, waveElevation, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('m');
title(['Wave Elevation for wave direction \beta_{wave} = ', ...
    num2str(rad2deg(beta_wave)), '°']);
grid on;

figure(2); clf;
% Plot the 6-DOF wave-frequency motions
DOF_txt = {'Surge (m)', 'Sway (m)', 'Heave (m)',...
    'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, T_scale(DOF) * eta_WF(:, DOF), 'LineWidth', 2);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF});
end

if ~isoctave
    sgtitle(['Wave-frequency motions for \beta_{wave} = ' ...
        num2str(rad2deg(beta_wave)), '°']);
end

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