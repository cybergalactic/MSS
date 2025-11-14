function [tau_wave1, waveElevation] = waveForceRAO(...
    t, S_M, Amp, Omega, mu, vessel, U, psi, beta_wave, numFreqIntervals)
% waveForceRAO computes the wave elevation and the 6-DOF generalized 1st-
% order wave forces, tau_wave1, on a ship using different wave spectra 
% (Modified Pierson-Moskowitz, JONSWAP, Torsethaugen) and Response 
% Amplitude Operators (RAOs) (Fossen 2021, Chapters 10.2.1 and 10.2.4). 
% The real and imaginary parts of the RAO tables are interpolated in 
% frequency and varying relative directions to compute the RAO amplitudes 
% and phases. This approach avoids unwrapping problems and interpolation 
% issues in RAO phase angles.
%
% INPUTS:
%   t                - Time vector (s)
%   S_M              - Spectral density matrix
%   Amp              - Wave amplitude matrix
%   Omega            - Wave frequencies (rad/s)
%   mu               - Wave spreading angles (radians)
%   vessel           - Vessel data structure containing RAO info
%   U                - Vessel speed (m/s)
%   psi              - Vessel heading angle (radians)
%   beta_wave        - Wave direction, 0 following sea, pi head sea (radians)
%   numFreqIntervals - Number of frequency intervals (> 50)
%
% OUTPUTS:
%   tau_wave1        - 6x1 generalized 1st-order wave forces (6-DOF)
%   waveElevation    - Wave elevation (m)
%
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-15
% Revisions:

persistent RAO_re_values_interpolated RAO_im_values_interpolated randomPhases;

% Constants
g = vessel.main.g; 

% Vesssel RAO data
freqs = vessel.forceRAO.w;   % RAO wave frequencies (rad/s)
raoAngles = deg2rad(vessel.headings); % RAO wave direction angles (rad)

% Initial computation and interpolation of real and imaginary RAO tables
if isempty(RAO_im_values_interpolated)

    numAngles = length(raoAngles);

    % Convert RAO amplitudes and phases to real and imaginary components
    RAO_re = cell(6, numAngles);
    RAO_im = cell(6, numAngles);
    for DOF = 1:6
        for j = 1:numAngles
            % Extract amp and phase for zero speed (index 1)
            % vessel.forceRAO.amp{DOF}(Omega, Dirctions, Speed)
            % vessel.forceRAO.phase{DOF}(Omega, Dirctions, Speed)
            RAO_amp = vessel.forceRAO.amp{DOF}(:,j,1);            
            RAO_phase = vessel.forceRAO.phase{DOF}(:,j,1);

            % Calculate the real and imaginary parts
            RAO_re{DOF, j} = RAO_amp .* cos(RAO_phase);
            RAO_im{DOF, j} = RAO_amp .* sin(RAO_phase);
        end
    end

    % Interpolate Re and Im parts of RAO to be valid for all Omega values.
    % Repeat this for all wave directions k = 1:numAngles
    rng(12345,"twister")
    numDirections = length(mu);
    randomPhases = 2 * pi * rand(numFreqIntervals, numDirections);

    RAO_re_values_interpolated = cell(1, 6);
    RAO_im_values_interpolated = cell(1, 6);
    for DOF = 1:6
        RAO_re_values = zeros(numFreqIntervals, numAngles);
        RAO_im_values = zeros(numFreqIntervals, numAngles);

        for k = 1:numAngles
            % Interpolate Re and Im parts of RAO for all Omega values
            RAO_re_values(:, k) = interp1(freqs, RAO_re{DOF, k}, Omega, ...
                'linear', 'extrap');
            RAO_im_values(:, k) = interp1(freqs, RAO_im{DOF, k}, Omega, ...
                'linear', 'extrap');
        end

        % Store the interpolated values for this DOF
        RAO_re_values_interpolated{DOF} = RAO_re_values;
        RAO_im_values_interpolated{DOF} = RAO_im_values;

    end

end

% Wave direction relative ship, beta_wave = 0 for following sea
beta_relative = beta_wave - psi;  

% Vector of spreading angles, scalar for M = 1 corresponding to mu = 0
beta_RAO = mod(beta_relative + mu, 2*pi); % Wrap to 0 to 2*pi

% Encounter frequency Omega_e(Omega, mu) for all frequencies and directions
Omega_e = abs(Omega - (Omega.^2 / g) * U .* cos(beta_RAO'));

% Compute the wave elevation (Fossen 2021, Eq. 10.83) using
% Amp = sqrt(2 * S_M * deltaOmega * deltaDirections).
% The summation over dim. 1 is frequencies and dim. 2 is directions 
if size(S_M, 2) == 1 
    % No spreading function, Amp(Omega) is a column vector
    waveElevation = sum( Amp .* cos(Omega_e * t + randomPhases), 1 );
else 
    % Directional wave spectrum, Amp(Omega, mu) is a matrix 
    waveElevation = sum( sum( Amp .* cos(Omega_e * t + randomPhases), 2 ), 1);
end

% Compute the complex RAOs as a function of frequency and wave direction
RAO_complex = cell(1, 6); % Initialize cell arrays
tau_wave1 = zeros(6,1);
numDirections = length(mu);
for DOF = 1:6

    % Retrieve stored frequency interpolated values
    RAO_re_values = RAO_re_values_interpolated{DOF};
    RAO_im_values = RAO_im_values_interpolated{DOF};

    % Initialize tables to stor the wave-direction interpolated results
    RAO_re_dir_interp = zeros(numFreqIntervals, numDirections);
    RAO_im_dir_interp = zeros(numFreqIntervals, numDirections);

    % Interpolate Re and Im parts of RAO for time-varying 'beta_RAO'
    % directions between 0 to 2*pi
    for k = 1:numDirections
        RAO_re_dir_interp(:, k) = interp1(raoAngles, RAO_re_values', ...
            beta_RAO(k), 'linear', 'extrap')';
        RAO_im_dir_interp(:, k) = interp1(raoAngles, RAO_im_values', ...
            beta_RAO(k), 'linear', 'extrap')';        
    end

    % Combine real and imaginary parts to form the complex RAO
    RAO_complex{DOF} = RAO_re_dir_interp + 1i * RAO_im_dir_interp;

    % Compute the generalized 1st-order wave forces (Fossen 2021, Eq. 10.96).
    if size(S_M, 2) == 1 
        % No spreading function/directional spectrum
        tau_wave1(DOF) = sum( abs(RAO_complex{DOF}) .* Amp .* ...
            cos(Omega_e * t + angle(RAO_complex{DOF}) + randomPhases), 1);
    else 
        % Directional spectrum
        tau_wave1(DOF) = sum( sum( abs(RAO_complex{DOF}) .* Amp .* ...
            cos(Omega_e * t + angle(RAO_complex{DOF}) + randomPhases), 2), 1);
    end
end

end


