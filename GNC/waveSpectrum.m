function [S, Omega, Amp] = waveSpectrum(spectrumType, spectrumParameters, ...
    numFreqIntervals, omegaMax)
% waveSpectrum is compatibel with MATLAB and GNU Octave (www.octave.org).
% Thos function generates a wave spectrum based on the specified type and 
% parameters.
% 
% [S, Omega, Amp] = waveSpectrum(spectrumType, spectrumParameters, ...
%     numFreqIntervals, omegaMax)
%
% INPUTS:
%   spectrumType       - Type of wave spectrum 
%                        ('Modified PM', 'JONSWAP', or 'Torsethaugen')
%   spectrumParameters - Parameters for the wave spectrum [Hs, w0] or [Hs, w0, gamma]
%                         Hs: Significant wave height
%                         w0: Peak frequency
%                         gamma: Peak enhancement factor (only for 'JONSWAP')
%   numFreqIntervals   - Number of frequency intervals (> 100)
%   omegaMax           - Maximum frequency (typicallly > 4 rad/s)
%
% OUTPUTS:
%   S                  - Wave spectrum (wave elevation power spectral 
%                        density function)
%   Omega              - Random frequency in each interval around the midpoint,
%                        confined to [-deltaOmega/2, deltaOmega/2]
%   Amp                - Wave amplitudes
%
% The function generates random frequencies within each interval deltaOmega
% and calculates the wave spectrum and amplitudes based on the given 
% spectrum type.
% 
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-06
% Revisions: 

Hs = spectrumParameters(1);
w0 = spectrumParameters(2);

if length(spectrumParameters) == 3
    gamma = spectrumParameters(3);
end

T0 = 2*pi / w0;
Tz = 0.710 * T0;

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

disp(['Significant wave height: Hs = ', num2str(Hs), ' m'])
disp(['Zero crossing period: Tz = ', num2str(Tz), ' s'])
disp(['Peak period: T0 = Tz / 0.710 = ', num2str(T0), ' s'])
disp(['Peak frequency: w0 = 2*pi/T0 = ', num2str(w0), ' rad/s'])
disp(['Number of wave frequensies: ', num2str(numFreqIntervals)])

if ~strcmp(spectrumType, 'Modified PM')
    % Modified PM spectrum
    S = wavespec(5, [Hs, Tz], Omega,0);
elseif ~strcmp(spectrumType, 'JONSWAP')
    % JONSWAP spectrum
    if w0 * sqrt(Hs) < 0 || w0 * sqrt(Hs) > 1.75
        error('It is reccomended to use 1.25 < (w0 * sqrt(Hs) < 1.75')
    end
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


end
