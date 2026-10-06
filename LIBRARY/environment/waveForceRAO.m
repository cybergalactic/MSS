function [tau_wave1, waveElevation] = waveForceRAO(...
    t, waveModel, U, psi, beta_wave)
% waveForceRAO computes the wave elevation and the 6-DOF generalized 1st-
% order wave forces, tau_wave1, on a ship using different wave spectra 
% (Modified Pierson-Moskowitz, JONSWAP, Torsethaugen) and Response 
% Amplitude Operators (RAOs) (Fossen 2027, Chapters 10.2.1 and 10.2.4).
% The real and imaginary parts of the RAO tables are interpolated in 
% frequency and varying relative directions to compute the RAO amplitudes 
% and phases. This approach avoids unwrapping problems and interpolation 
% issues in RAO phase angles.
%
% INPUTS:
%   t          - Time (s)
%   waveModel  - Wave model returned by waveInitialization
%   U          - Vessel speed (m/s)
%   psi        - Vessel heading angle (rad)
%   beta_wave  - Wave direction, 0 following sea, pi head sea (rad)
%
% OUTPUTS:
%   tau_wave1        - 6x1 generalized 1st-order wave forces (6-DOF)
%   waveElevation    - Wave elevation (m)
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-15
% Revisions:
%   2026-10-03: Added periodic RAO interpolation at 0/2*pi radians.
%   2026-10-03: Use the caller-controlled random-number stream.
%   2026-10-03: Replaced persistent state with an explicit waveModel.

if ~strcmp(waveModel.raoType, 'force')
    error('waveForceRAO:InvalidWaveModel', ...
        'waveModel must be initialized with force RAOs.');
end

S_M = waveModel.S_M;
Amp = waveModel.Amp;
Omega = waveModel.Omega;
mu = waveModel.mu;
g = waveModel.g;
randomPhases = waveModel.randomPhases;
numFreqIntervals = length(Omega);

% Wave direction relative ship, beta_wave = 0 for following sea
beta_relative = beta_wave - psi;  

% Vector of spreading angles, scalar for M = 1 corresponding to mu = 0
beta_RAO = mod(beta_relative + mu, 2*pi); % Wrap to 0 to 2*pi

% The initialization closes the directional RAO tables at 2*pi.
raoAngles = waveModel.raoAngles;

% Encounter frequency Omega_e(Omega, mu) for all frequencies and directions
Omega_e = abs(Omega - (Omega.^2 / g) * U .* cos(beta_RAO'));

% Compute the wave elevation (Fossen 2027, Eq. 10.83) using
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

    % Retrieve the initialized, periodically closed force-RAO tables.
    RAO_re_values = waveModel.RAO_re{DOF};
    RAO_im_values = waveModel.RAO_im{DOF};

    % Initialize tables to store the wave-direction interpolated results
    RAO_re_dir_interp = zeros(numFreqIntervals, numDirections);
    RAO_im_dir_interp = zeros(numFreqIntervals, numDirections);

    % Interpolate Re and Im parts of RAO for time-varying 'beta_RAO'
    % directions between 0 to 2*pi
    for k = 1:numDirections
        RAO_re_dir_interp(:, k) = interp1(raoAngles, ...
            RAO_re_values', beta_RAO(k), 'linear')';
        RAO_im_dir_interp(:, k) = interp1(raoAngles, ...
            RAO_im_values', beta_RAO(k), 'linear')';
    end

    % Combine real and imaginary parts to form the complex RAO
    RAO_complex{DOF} = RAO_re_dir_interp + 1i * RAO_im_dir_interp;

    % Compute the generalized 1st-order wave forces (Fossen 2027, Eq. 10.96).
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


