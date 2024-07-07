function [S_M, Omega, Amp, S, M, mu] = waveSpectrum(spectrumType, ...
    Parameters, numFreqIntervals, omegaMax, spreadingFlag, numDirections)
% waveSpectrum is compatible with MATLAB and GNU Octave (www.octave.org). 
% This function generates a directional wave spectrum 
% 
%   S_M(Omega, mu) = S(Omega) * M(mu)
% 
% where 'Omega' is the frequency vector and 'mu' is wave direction vector.
% The wave amplitudes 'Amp' for different spectrum types 'S(Omega)' and 
% directional wave spectrum 'M(mu)' are also computed. 
% 
% The optional flag 'spreadingFlag' allows for the propagation of waves in 
% multiple directions around a main propagation direction (mu = 0). When 
% enabled, it spreads the energy over 'numDirections' wave directions 
% within an angle range of [-pi/2, pi/2] from the wind direction, providing 
% a more realistic representation of natural sea states.
% 
% [S, Omega, Amp, M, mu] = waveSpectrum(spectrumType, Parameters, ...
%    numFreqIntervals, omegaMax, spreadingFlag, numDirections)
%
% INPUTS:
%   spectrumType     - Type of wave spectrum, S(Omega),
%                      ('Modified PM', 'JONSWAP', or 'Torsethaugen')
%   Parameters       - Parameters for the wave spectrum [Hs, w0] or 
%                      [Hs, w0, gamma], where
%                         Hs: Significant wave height
%                         w0: Peak frequency
%                         gamma: Peak enhancement factor (only for 'JONSWAP')
%   numFreqIntervals - Number of frequency intervals (> 100)
%   omegaMax         - Maximum frequency (typicallly > 4 rad/s)
%   spreadingFlag    - Optional flag for wave spreading (default is false)
%   numDirections    - Number of wave directions (default is 36)
%
% OUTPUTS:
%
%   S_M              - Matrix of wave spectrum with spreading 
%                      (frequencies x directions)
%   Omega            - Random frequency in each interval around the midpoint,
%                      confined to [-deltaOmega/2, deltaOmega/2]
%   Amp              - Wave amplitudes
%   S                - Optional wave spectrum S(Omega), wave elevation 
%                      power spectral density function
%   M                - Optional wave directional spectrum, M(mu)
%   mu               - Optional wave directions
%
% The function generates random frequencies within each interval deltaOmega
% and calculates the wave spectrum and amplitudes based on the given 
% spectrum type.
%
% Examples:
%   [S_M, Omega] = waveSpectrum('Modified PM', [Hs, w0])
%   [S_M, Omega] = waveSpectrum('JONSWAP', [Hs, w0, gamma])
%   [S_M, Omega] = waveSpectrum('Torsethaugen', [Hs, w0])
%   [S_M, Omega, Amp] = waveSpectrum('Modified PM', [Hs, w0], true)
%   [S_M, Omega, Amp] = waveSpectrum('Modified PM', [Hs, w0], true, 24)
% 
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-07-06
% Revisions: 

if nargin < 6
    numDirections = 16; % Default value for number of wave directions
end

if nargin < 5
    spreadingFlag = false;                                                                                                                        false; % Default value for spreadingFlag
end

Hs = Parameters(1);
w0 = Parameters(2);

if length(Parameters) == 3
    gamma = Parameters(3);
end

if spreadingFlag
    % Calculate the spreading function M(mu) and normalize it
    mu = linspace(-pi/2, pi/2, numDirections)'; % Wave directions
    M = (2/pi) * cos(mu).^2; % Spreading function, mu in [-pi/2, pi/2]
    integral_M = trapz(mu, M); % Integral of M using trapezoidal rule
    M = M / integral_M; % Normalize to ensure the integral is 1
else
    % No spreading
    M = 1; 
    mu = 0;  
    numDirections = 1;
end

T0 = 2*pi / w0;
Tz = 0.710 * T0;

% Generate fixed frequencies within the range [0, omegaMax]
Omega = linspace(eps, omegaMax, numFreqIntervals)';

if M == 1
    disp(['Spectrum: ', spectrumType, ' (no spreading function)'])
else
    disp(['Spectrum: ', spectrumType, ' with spreading function M(mu)'])
end

if strcmp(spectrumType, 'Modified PM')
    % Modified PM spectrum
    S = wavespec(5, [Hs, Tz], Omega,0);
elseif strcmp(spectrumType, 'JONSWAP')
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

disp(['Significant wave height: Hs = ', num2str(Hs), ' m'])
disp(['Zero crossing period: Tz = ', num2str(Tz), ' s'])
disp(['Peak period: T0 = Tz / 0.710 = ', num2str(T0), ' s'])
disp(['Peak frequency: w0 = 2*pi/T0 = ', num2str(w0), ' rad/s'])
disp(['Number of wave frequensies: ', num2str(numFreqIntervals)])
disp(['Number of wave directions: ', num2str(numDirections)])

% Compute the 2-D wave spectrum S_M(Omega, mu) = S(Omega) * M(mu)
S_M = S(:) * M(:).';

% Compute the wave amplitudes
Amp = sqrt(2 * S_M);

end
