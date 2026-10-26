function [componentPhase, waveModel] = waveEncounterPhase(t, waveModel, Omega_e)
% waveEncounterPhase accumulates signed encounter phase (Fossen 2027, Ch. 10).
% Omega_e is the signed encounter rate for each frequency-direction pair:
%   Omega_e = Omega - (Omega.^2/g)*U*cos(beta_wave - psi + mu).
% The nonnegative encounter-frequency magnitude is abs(Omega_e), but phase
% is advanced using the signed rate, including overtaking in following seas.
%
% Return waveModel to the next RAO call to preserve the phase history.
% Successive rates are integrated by the trapezoidal rule. On the first
% call, the current rate is assumed constant from time zero to t. Calling
% again at the same time does not advance phase. To restart a realization,
% obtain a new waveModel from waveInitialization; no persistent state is used.
%
% Author: Thor I. Fossen
% Date:   2026-10-26

if ~isnumeric(t) || ~isreal(t) || ~isscalar(t) || ~isfinite(t) || t < 0
    error('waveEncounterPhase:InvalidTime', ...
        'Time must be a finite, nonnegative scalar.');
end

if ~isfield(waveModel, 'phaseTime') || isempty(waveModel.phaseTime)
    componentPhase = waveModel.randomPhases + t * Omega_e;
else
    dt = t - waveModel.phaseTime;
    if dt < 0
        error('waveEncounterPhase:TimeOrder', ...
            'Time must not decrease. Reinitialize waveModel to restart.');
    end
    componentPhase = waveModel.componentPhase + ...
        (dt/2) * (waveModel.encounterRate + Omega_e);
end

waveModel.componentPhase = componentPhase;
waveModel.encounterRate = Omega_e;
waveModel.phaseTime = t;
end
