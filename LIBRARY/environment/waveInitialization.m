function [environment, waveModel] = waveInitialization(vessel, environment, raoType)
% waveInitialization initializes a directional wave spectrum and RAO model.
% The function prepares the time-invariant data and initializes the phase
% history used by waveForceRAO or waveMotionRAO. Call it once before the
% simulation loop and pass waveModel to the selected RAO function.
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
%                 The frequency band is chosen automatically from the spectrum
%                 and the available RAO data; numFreqIntervals is unchanged.
%   raoType      - Optional RAO model: 'force' (default) or 'motion'
%
% Outputs:
%   environment - Input structure augmented with S_M, Amp, Omega, mu,
%                 frequencyBand, and retainedEnergyFraction
%   waveModel   - Initialized wave spectrum, random phases, and interpolated
%                 RAO tables used by waveForceRAO or waveMotionRAO. The RAO
%                 functions return an updated waveModel containing encounter
%                 phase history; pass it to the next call. For a
%                 submerged vessel, raoReferenceDepth is the depth at which
%                 the hydrodynamic data were generated.
%
% The random-number generator is intentionally not set here. Use rng at the
% simulation or example level before calling waveInitialization when a
% repeatable wave realization is required.
% The automatic band targets the central 99% of spectral energy, clipped to
% the finite RAO range. A warning reports retention below 95% or a frequency
% interval containing over 10% of retained energy. Amplitudes are not rescaled
% to compensate for truncation. The component count remains user-selected.
%
% Author:    Thor I. Fossen
% Date:      2026-10-03
% Revisions:
%   2026-10-14: Store the RAO reference depth for submerged-vessel loads.
%   2026-10-26: Initialize explicit encounter-phase history for maneuvering.
%   2026-10-27: Select an energy-based frequency band within the RAO data.

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

% Use the finite RAO grid. ShipX/VERES uses 10 rad/s as an infinity marker;
% that convention must not remove genuine high-frequency data from other solvers.
allFreqs = raoData.w(:);
freqIndex = isfinite(allFreqs) & allFreqs >= 0;
if isfield(vessel,'roll') && isfield(vessel.roll,'veres')
    freqIndex = freqIndex & abs(allFreqs - 10) > 10*eps(10);
end
freqs = allFreqs(freqIndex);
if numel(freqs) < 2 || any(diff(freqs) <= 0)
    error('waveInitialization:InvalidFrequencies', ...
        'At least two strictly increasing finite RAO frequencies are required.');
end
validateattributes(environment.w0, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(environment.Hs, {'numeric'}, {'scalar','real','finite','nonnegative'});
if environment.w0 <= freqs(1) || environment.w0 >= freqs(end)
    error('waveInitialization:PeakOutsideRAORange', ...
        'Peak frequency %.3g rad/s must lie inside the RAO range [%.3g, %.3g] rad/s.', ...
        environment.w0, freqs(1), freqs(end));
end

[frequencyBand, retainedEnergy] = spectrumBand(environment, freqs);
environment.frequencyBand = frequencyBand;
environment.retainedEnergyFraction = retainedEnergy;
if isfield(environment,'maxFreq')
    environment = rmfield(environment,'maxFreq'); % Retire the manual cutoff.
end

% Directional spectrum and component amplitudes
[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum( ...
    environment.spectrumType, spectrumParameters, ...
    environment.numFreqIntervals, frequencyBand(2), ...
    environment.spreadingFlag, environment.numDirections, frequencyBand(1));

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
waveModel.raoReferenceDepth = 0;
if isfield(vessel.main, 'submerged') && ...
        isscalar(vessel.main.submerged) && logical(vessel.main.submerged) && ...
        isfield(vessel.main, 'submergenceDepth') && ...
        isscalar(vessel.main.submergenceDepth) && ...
        isfinite(vessel.main.submergenceDepth)
    waveModel.raoReferenceDepth = max(0, vessel.main.submergenceDepth);
end

numFreqIntervals = length(Omega);
numDirections = length(mu);
waveModel.randomPhases = 2 * pi * rand(numFreqIntervals, numDirections);
waveModel.componentPhase = waveModel.randomPhases;
waveModel.encounterRate = [];
waveModel.phaseTime = [];

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
            'linear');
        RAO_im_values(:, k) = interp1(freqs, RAO_im, Omega, ...
            'linear');
    end

    % Close each table periodically: the RAO at 2*pi equals the RAO at 0.
    waveModel.RAO_re{DOF} = [RAO_re_values, RAO_re_values(:, 1)];
    waveModel.RAO_im{DOF} = [RAO_im_values, RAO_im_values(:, 1)];
end

end

function [band, retained] = spectrumBand(environment, freqs)
% Integrate on a dense reference grid, independently of the user's component
% count. The central 99% energy band avoids spending components on weak tails.
if environment.Hs == 0
    % With no wave energy, spectral quantiles are undefined.
    band = [freqs(1),freqs(end)];
    retained = 1;
    return
end
w0 = environment.w0;
w = logspace(log10(w0/100),log10(100*w0),16385)';
switch environment.spectrumType
    case 'Modified PM'
        S = waveSpectrum(5,[1,0.710*2*pi/w0],w,0);
    case 'JONSWAP'
        S = waveSpectrum(7,[1,w0,3.3],w,0);
    case 'Torsethaugen'
        % Its shape depends on Hs as well as w0.
        S = waveSpectrum(8,[environment.Hs,w0],w,0);
    otherwise
        error('waveInitialization:UnknownSpectrum','Unknown wave spectrum type.');
end
energy = cumtrapz(w,S(:));
if any(~isfinite(energy)) || energy(end) <= 0
    error('waveInitialization:InvalidSpectrum','The spectrum must have finite positive energy.');
end
cdf = energy / energy(end);
[cdfUnique, indices] = unique(cdf);
band = interp1(cdfUnique,w(indices),[0.005 0.995]);
band = [max(band(1),freqs(1)), min(band(2),freqs(end))];
if band(1) >= band(2)
    error('waveInitialization:NoSpectralOverlap', ...
        'The energetic spectrum band does not overlap the available RAO range.');
end
retained = diff(interp1(w,cdf,band));
if retained < 0.95 && environment.Hs > 0
    warning('waveInitialization:SpectrumTruncated', ...
        ['RAO range limits the wave band to [%.3g, %.3g] rad/s, retaining %.1f%% ' ...
         'of spectral energy. The requested Hs is not fully represented; ' ...
         'the spectrum is not rescaled.'],band(1),band(2),100*retained);
end
% Check frequency quadrature without changing the requested component count.
n = environment.numFreqIntervals;
validateattributes(n, {'numeric'}, {'scalar','integer','positive','finite'});
edges = linspace(band(1),band(2),n+1);
fractions = diff(interp1(w,cdf,edges));
if max(fractions) > 0.1*retained && environment.Hs > 0
    warning('waveInitialization:SpectrumUnderresolved', ...
        ['One frequency interval contains over 10%% of retained spectral energy. ' ...
         'Consider increasing numFreqIntervals (currently %d).'],n);
end
end
