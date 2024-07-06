function exWaveFrequencyMotion()
% exWaveFrequencyMotion is compatibel with MATLAB and GNU Octave (www.octave.org).
% This function computes the wave elevation and the wave-frequency (WF) 
% motion, eta_w, on a ship using different wave spectra (Modified Pierson-
% Moskowitz, JONSWAP, Torsethaugen) and Response Amplitude Operator (RAO) 
% data (Fossen 2021, Chapters 10.2.1 and 10.2.3). The total motion is:
%
%   y = eta + eta_w
%
% where eta is the low-frequency ship model component. The real and 
% imaginary parts of the RAO tables are interpolated in frequency for 
% varying wave directions to compute the RAO amplitudes and phases. This 
% approach avoids unwrapping problems and interpolation issues in RAO 
% phase angles.
displayMfileHeader('exWaveFrequencyMotion.m');

% Call the function to select ship and load mat-file containing the motion RAOs
matFile = selectShip(); 
load(which(matFile), 'vessel');
disp(['Loaded the motion RAO structure, vessel.motionRAO, from ', matFile])

% Call the function to select wave spectrum
[noSpectrum, nameSpectrum] = selectWaveSpectrum(); 

%% Simulation parameters
Tfinal = 100;                   % Duration of the simulation (s)
h = 0.05;                       % Time step (s)
numFreqIntervals = 100;         % Number of frequency intervals (>100)
g = 9.81;                       % Acceleration of gravity (m/s^2)

% Sea state - JONSWAP spectrum with DNV (2007) correction
Hs = 10;                        % Significant wave height (m)
Tz = 10;                        % Zero-crossing period (s)
    
% Wave direction relative bow, 0 deg for following sea, 180 deg for head sea
beta_wave = deg2rad(140);

% Ship (can be time varying)
U = 0;                          % Speed (m/s)
psi = deg2rad(0);               % Heading angle (rad)

% Fixed seed for reproducibility
seed = 12345; 
rng(seed);

freqs = vessel.motionRAO.w;             % RAO wave frequencies
raoAngles = vessel.headings;            % RAO wave direction angles
numAngles = length(raoAngles);          
omegaMax = freqs(end);

% Random phases (Fossen 2021, Eq. 10.76)
randomPhases = 2 * pi * rand(numFreqIntervals, 1);

% Random frequencies within each interval deltaOmega (Fossen 2021, Eq. 10.77)
deltaOmega = omegaMax / numFreqIntervals; % Frequency interval
Omega = zeros(numFreqIntervals, 1);
for i = 1:numFreqIntervals
    % Midpoint of the current interval
    omega_mid = (i-1) * deltaOmega + deltaOmega / 2;
    % Random frequency in each interval around the midpoint,
    % confined to [-deltaOmega/2, deltaOmega/2]
    Omega(i) = omega_mid + (rand() - 0.5) * deltaOmega;
end

% Calculate the wave spectrum power intensity S(Omega) for each frequency
T0 = Tz / 0.710;
w0 = 2*pi / T0;  % Modal (peak) frequency

disp(['Number of wave frequensies: ', num2str(numFreqIntervals)])
disp(['Wave direction: ', num2str(rad2deg(beta_wave)), '°'])
disp(' ');
disp(['Significant wave height: Hs = ', num2str(Hs), ' m'])
disp(['Zero crossing period: Tz = ', num2str(Tz), ' s'])
disp(['Peak period: T0 = Tz / 0.710 = ', num2str(T0), ' s'])
disp(['Peak frequency: w0 = 2*pi/T0 = ', num2str(w0), ' rad/s'])

if noSpectrum == 1
    % Modified PM spectrum
    S = wavespec(5, [Hs, Tz], Omega,0); 
elseif noSpectrum  == 2
    % JONSWAP spectrum
    if w0 * sqrt(Hs) < 0 || w0 * sqrt(Hs) > 1.75
        error('It is reccomended to use 1.25 < (w0 * sqrt(Hs) < 1.75')
    end
    gamma = 3.3;
    S = wavespec(7, [Hs, w0, gamma], Omega, 0); 
else
    % Torsethaugen spectrum
    if w0 < 0.6
        disp('For w0 <= 0.6, only one peak in the Torsetaugen spectrum appears')
    end
    S = wavespec(8, [Hs, w0], Omega,0); 
end

% Calculate the wave amplitudes (Fossen 2021, Eq. 10.75)
Amp = sqrt(2 * S * deltaOmega);

% Convert RAO amplitudes and phases to real and imaginary components
RAO_re = cell(6, numAngles);
RAO_im = cell(6, numAngles);
for k = 1:6
    for j = 1:numAngles
        % Extract amp and phase for zero speed, index 1
        RAO_phase = vessel.motionRAO.phase{k}(:,:,1);
        RAO_amp = vessel.motionRAO.amp{k}(:,:,1);

        % Calculate the real and imaginary parts
        RAO_re{k, j} = RAO_amp .* cos(RAO_phase);
        RAO_im{k, j} = RAO_amp .* sin(RAO_phase);
    end
end

% Compute the relative wave direction for the current heading
beta_relative = beta_wave - psi;
beta_relative = mod(beta_relative, 2*pi); % Wrap to 0 and 2*pi

% Compute the encounter frequency Omega_e (Fossen 2021, Eq. 10.85)
Omega_e = abs(Omega - (Omega.^2 / g) * U * cos(beta_relative));

%% Interpolate RAOs in real-time
RAO_complex = cell(1, 6); % Initialize cell arrays
for DOF = 1:6
    % Interpolate RAOs for each frequency
    F_re_values = zeros(numAngles, numFreqIntervals);
    F_im_values = zeros(numAngles, numFreqIntervals);

    for k = 1:numAngles
        % Interpolate Re and Im parts of RAO for frequencies Omega
        F_re_values(k, :) = interp1(freqs(:), RAO_re{DOF, k}(:,k), ...
            Omega, 'linear', 'extrap');
        F_im_values(k, :) = interp1(freqs(:), RAO_im{DOF, k}(:,k), ...
            Omega, 'linear', 'extrap');
    end

    % Interpolate Re and Im parts of RAO for wave directions 0 to 2*pi
    F_re_dir_interp = interp1(raoAngles, F_re_values, beta_relative,...
        'linear', 'extrap');
    F_im_dir_interp = interp1(raoAngles, F_im_values, beta_relative,...
        'linear', 'extrap');

    % Combine real and imaginary parts to form the complex RAO
    RAO_complex{DOF} = F_re_dir_interp + 1i * F_im_dir_interp;
end

%% MAIN LOOP
t = 0:h:Tfinal; % Time vector
eta_WF = zeros(length(t), 6); 
waveElevation = zeros(length(t), 1);

% Compute the wave elevation and the 6-DOF wave-frequency (WF) motions
for t_i = 1:length(t)
    % Wave elevation (Fossen 2021, Eq. 10.76)
    waveElevation(t_i) = sum(Amp .* cos(Omega_e .* t(t_i) + randomPhases), 1);

    % 6-DOF wave-frequency motions (Fossen 2021, Eq. 10.105)
    for DOF = 1:6
        eta_WF(t_i, DOF) = sum(abs(RAO_complex{DOF})' .* Amp ...
            .* cos(Omega_e .* t(t_i) + angle(RAO_complex{DOF})' + randomPhases));
    end
end


%% PLOTS
figure(1);
% Plot the wave spectrum
subplot(311);
plot(Omega, S, 'b-', 'LineWidth', 2);
xlabel('Omega (rad/s)');
ylabel('m^2 s');
title(nameSpectrum);
xlim([0, omegaMax]);
hold on; 
plot([w0, w0], [min(S), max(S)], 'r', 'LineWidth', 2);
hold off; 
grid on;

% Plot the wave amplitude
subplot(312);
plot(Omega, Amp, 'b-', 'LineWidth', 2);
xlabel('Omega (rad/s)');
ylabel('m');
title('Wave Amplitudes');
xlim([0, omegaMax]);
grid on;

% Plot the wave elevation
subplot(313);
plot(t, waveElevation, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('m');
title(['Wave Elevation for wave direction \beta_{wave} = ', ...
    num2str(rad2deg(beta_wave)), '°']);
grid on;

figure(2);
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
function matFile = selectShip()

% Define ship options
ships = {'Supply vessel', 'S175 container ship', 'Tanker', ...
    'FPSO', 'Semisubmersible'};

% Display menu to user and get selection
fprintf('Please choose a ship:\n');
for i = 1:length(ships)
    fprintf('  %d. %s\n', i, ships{i});
end
choice = input('  Enter the number corresponding to your choice: ');

% Validate user input
while choice < 1 || choice > length(ships)
    disp('  Invalid choice. Please try again.');
    choice = input('  Enter the number corresponding to your choice: ');
end

switch choice
    case 1
        matFile = 'supply.mat';
    case 2
        matFile = 's175.mat';
    case 3
        matFile = 'tanker.mat';
    case 4
        matFile = 'FPSO.mat';
    case 5
        matFile = 'semisub.mat';
    otherwise
        disp('  Error: Unknown ship type.');
end

disp(' ');

end

% Function to select wave spectrum
function [noSpectrum, nameSpectrum] = selectWaveSpectrum()
    spectra = {'Modifed Pierson-Moskowitz (PM) spectrum', ...
        'JONSWAP spectrum', 'Torsethaugen spectrum'};
    
    % Display menu to user and get selection
    fprintf('\nPlease choose a wave spectrum:\n');
    for i = 1:length(spectra)
        fprintf('  %d. %s\n', i, spectra{i});
    end
    choice = input('  Enter the number corresponding to your choice: ');

    % Validate user input
    while choice < 1 || choice > length(spectra)
        disp('  Invalid choice. Please try again.');
        choice = input('  Enter the number corresponding to your choice: ');
    end
    
    noSpectrum = choice;
    nameSpectrum = spectra{noSpectrum};

    disp(' ');
    disp(nameSpectrum);

end

