function [tau_wave2, waveModel] = waveDriftRAO(waveModel, psi, beta_wave, depth)
% waveDriftRAO computes 6-DOF generalized second-order mean wave-drift loads.
% The calculation uses the real, signed diagonal drift coefficients prepared
% by waveInitialization and the squared component amplitudes of a long- or
% short-crested irregular sea (Fossen 2027, Chapters 5 and 10):
%
%   tau_wave2(dof) = sum H_wave2(dof,Omega,beta_RAO) .* Amp.^2
%
% The MSS vessel.drift layout contains coefficients for physical DOFs 1, 2,
% and 6 (surge, sway, and yaw). Consequently, this function returns zero in
% heave, roll, and pitch. It does not infer a load application point. If the
% dynamics origin differs from the hydrodynamic-data origin, transform the
% complete generalized load explicitly, for example
%
%   tau_CO = Hmtrx(r_CO_to_data)' * tau_wave2
%
% Such a reference-point transformation can create roll and pitch moments
% from the horizontal drift forces. Any stored yaw coefficient is retained.
%
% INPUTS:
%   waveModel - Force- or motion-RAO model returned by waveInitialization.
%               It must contain vessel.drift data prepared at initialization.
%   psi       - Craft heading angle (rad)
%   beta_wave - Wave propagation direction in NED (rad). Relative direction
%               beta_wave - psi is 0 in following seas, pi in head seas.
%   depth     - Optional submerged-body depth, positive down (m). When
%               supplied, each squared-amplitude drift contribution is
%               attenuated below the RAO reference depth by exp(-2*k*Delta_d), 
%               where k = Omega^2/g. Omit this input for the original 
%               surface-vessel behavior.
%
% OUTPUTS:
%   tau_wave2 - 6x1 generalized mean wave-drift load at the hydrodynamic-data
%               reference origin: [X, Y, 0, 0, 0, N]' (N, N, N, Nm, Nm, Nm)
%   waveModel - Unchanged input model, returned for paired calling syntax with
%               waveForceRAO. Mean drift has no encounter-phase state.
%
% The coefficients are loads per squared incident-wave amplitude (N/m^2 for
% surge and sway, Nm/m^2 for yaw). Only the supplied zero-speed coefficient
% table is used. Forward speed affects neither these coefficients nor their
% heading argument; use speed-dependent drift data only after extending the
% data model explicitly. Equal-frequency mean drift is included, while
% fluctuating difference-frequency and sum-frequency QTF loads are not.
%
% Example:
%   [environment, waveModel] = waveInitialization(vessel, environment);
%   [tau_wave1, zeta, waveModel] = waveForceRAO( ...
%       t, waveModel, U, psi, beta_wave);
%   [tau_wave2, waveModel] = waveDriftRAO(waveModel, psi, beta_wave);
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-11-14
% Revisions:
%   2026-11-14: Add optional squared-amplitude depth attenuation for
%               submerged-body mean-drift loads.

if ~isstruct(waveModel) || ~isfield(waveModel, 'hasDriftRAO') || ...
        ~isscalar(waveModel.hasDriftRAO) || ~waveModel.hasDriftRAO || ...
        ~isfield(waveModel, 'driftCoeff')
    error('waveDriftRAO:InvalidWaveModel', ...
        ['waveModel must be returned by waveInitialization for a vessel ' ...
         'containing signed mean wave-drift coefficients.']);
end
validateattributes(psi, {'numeric'}, {'scalar','real','finite'}, ...
    mfilename, 'psi');
validateattributes(beta_wave, {'numeric'}, {'scalar','real','finite'}, ...
    mfilename, 'beta_wave');

Amp = waveModel.Amp;
mu = waveModel.mu;
raoAngles = waveModel.raoAngles;

% Mean drift scales with incident-wave amplitude squared. Apply the square
% of the first-order deep-water amplitude attenuation when depth is supplied.
depthAttenuation = 1;
if nargin >= 4 && ~isempty(depth)
    if ~isscalar(depth) || ~isnumeric(depth) || ~isreal(depth) || ...
            ~isfinite(depth)
        error('waveDriftRAO:InvalidDepth', ...
            'depth must be a finite real scalar in meters, positive down.');
    end
    referenceDepth = 0;
    if isfield(waveModel, 'raoReferenceDepth')
        referenceDepth = waveModel.raoReferenceDepth;
    end
    depthIncrement = max(depth - referenceDepth, 0);
    waveNumber = waveModel.Omega.^2 / waveModel.g;
    depthAttenuation = exp(-2 * waveNumber * depthIncrement);
end

% Use the same propagation-direction convention and periodic heading
% interpolation as waveForceRAO.
beta_relative = beta_wave - psi;
beta_RAO = mod(beta_relative + mu, 2*pi);

tau_wave2 = zeros(6, 1);
physicalDOFs = [1 2 6];
for j = 1:3
    coeffValues = waveModel.driftCoeff{j};
    coeffDirInterp = interp1(raoAngles, coeffValues', ...
        beta_RAO, 'linear')';
    components = depthAttenuation .* coeffDirInterp .* Amp.^2;
    tau_wave2(physicalDOFs(j)) = sum(components(:));
end

end
