function [environment, waveModel] = waveInitialization(vessel, environment, raoType)
% waveInitialization initializes a directional wave spectrum and RAO model.
% The function prepares all time-invariant data used by waveForceRAO or
% waveMotionRAO. Call it once before the simulation loop and pass waveModel
% to the selected RAO function.
%
% Inputs:
%   vessel      - Vessel structure containing main.g, headings, and forceRAO
%                 or motionRAO
%   environment - Wave-environment structure with the fields:
%                   Hs               significant wave height (m)
%                   w0               spectral peak frequency (rad/s)
%                   spectrumType     'Modified PM', 'JONSWAP', or
%                                    'Torsethaugen'
%                   spreadingFlag    true for a directional spectrum
%                   numFreqIntervals number of frequency intervals
%                   numDirections    number of spreading directions
%                 Optional field:
%                   maxFreq          maximum RAO frequency (rad/s),
%                                    default 5 rad/s
%   raoType      - Optional RAO model: 'force' (default) or 'motion'
%
% Outputs:
%   environment - Input structure augmented with S_M, Amp, Omega, and mu
%   waveModel   - Initialized wave spectrum, random phases, and interpolated
%                 RAO tables used by waveForceRAO or waveMotionRAO
%
% The random-number generator is intentionally not set here. Use rng at the
% simulation or example level before calling waveInitialization when a
% repeatable wave realization is required.
%
% Author:    Thor I. Fossen
% Date:      2026-10-03
% Revisions:

if nargin < 3
    raoType = 'force';
end

switch lower(raoType)
    case {'force', 'forcerao'}
        raoData = vessel.forceRAO;
        waveModel.raoType = 'force';
    case {'motion', 'motionrao'}
        raoData = vessel.motionRAO;
        waveModel.raoType = 'motion';
    otherwise
        error('waveInitialization:InvalidRAOType', ...
            'raoType must be ''force'' or ''motion''.');
end

% Wave-spectrum parameters
spectrumParameters = [environment.Hs, environment.w0];
if strcmp(environment.spectrumType, 'JONSWAP')
    spectrumParameters = [spectrumParameters, 3.3];
end

% Limit the spectrum to finite frequencies represented by the selected RAOs.
% A value of 10 rad/s in hydrodynamic data is commonly used for infinity.
requestedMaxFreq = 5;
if isfield(environment, 'maxFreq')
    requestedMaxFreq = environment.maxFreq;
end

allFreqs = raoData.w(:);
finiteFreqs = allFreqs(allFreqs < 10);
if isempty(finiteFreqs)
    error('waveInitialization:NoFiniteFrequencies', ...
        'The selected RAO model contains no finite frequencies below 10 rad/s.');
end

maxFreq = min(requestedMaxFreq, max(finiteFreqs));
freqIndex = allFreqs <= maxFreq;
freqs = allFreqs(freqIndex);
if length(freqs) < 2
    error('waveInitialization:InsufficientFrequencies', ...
        'At least two RAO frequencies are required.');
end

% Directional spectrum and component amplitudes
[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum( ...
    environment.spectrumType, spectrumParameters, ...
    environment.numFreqIntervals, freqs(end), ...
    environment.spreadingFlag, environment.numDirections);

environment.S_M = S_M;
environment.Omega = Omega;
environment.mu = mu;
environment.Amp = Amp;

% Store the wave quantities evaluated by the selected RAO function.
waveModel.S_M = S_M;
waveModel.Amp = Amp;
waveModel.Omega = Omega;
waveModel.mu = mu;
waveModel.g = vessel.main.g;

numFreqIntervals = length(Omega);
numDirections = length(mu);
waveModel.randomPhases = 2 * pi * rand(numFreqIntervals, numDirections);

% Convert the RAOs to real and imaginary parts and interpolate them once on
% the wave-frequency grid. Index 1 is the zero-speed RAO table.
raoAngles = vessel.headings(:);
numAngles = length(raoAngles);
waveModel.raoAngles = [raoAngles; 2*pi];
waveModel.RAO_re = cell(1, 6);
waveModel.RAO_im = cell(1, 6);

for DOF = 1:6
    RAO_re_values = zeros(numFreqIntervals, numAngles);
    RAO_im_values = zeros(numFreqIntervals, numAngles);

    for k = 1:numAngles
        RAO_amp = raoData.amp{DOF}(freqIndex, k, 1);
        RAO_phase = raoData.phase{DOF}(freqIndex, k, 1);
        RAO_re = RAO_amp .* cos(RAO_phase);
        RAO_im = RAO_amp .* sin(RAO_phase);

        RAO_re_values(:, k) = interp1(freqs, RAO_re, Omega, ...
            'linear', 'extrap');
        RAO_im_values(:, k) = interp1(freqs, RAO_im, Omega, ...
            'linear', 'extrap');
    end

    % Close each table periodically: the RAO at 2*pi equals the RAO at 0.
    waveModel.RAO_re{DOF} = [RAO_re_values, RAO_re_values(:, 1)];
    waveModel.RAO_im{DOF} = [RAO_im_values, RAO_im_values(:, 1)];
end

end
