function cfg = hydroVesselConfig(matFile, useGenericDefaults)
% cfg = hydroVesselConfig(matFile,useGenericDefaults) returns draft
% simulation, environment,
% controller, and damping parameters for the hydrodynamic vessel data file
% matFile. Common defaults are defined first and vessel-specific values are
% applied afterwards. If a vessel case does not explicitly override its
% linear damping parameters, they are loaded from vessel.powerBased in the
% MAT-file.
%
% Input:
%   matFile: Vessel data filename, for example 'supply.mat'
%   useGenericDefaults: Optional logical flag for a custom MAT-file. When
%                      true, editable conservative defaults are used instead
%                      of a named vessel preset.
%
% Output:
%   cfg: Configuration structure used by guiSIMhydroVessel and SIMhydroVessel
%
% Author:    Thor I. Fossen
% Date:      2026-10-04
% Revisions:
%   2026-10-04 Made heading autopilot the default except for the supply
%              vessel and semisubmersible, which use DP control.
%   2026-10-05 Use vessel.powerBased as the default source for Capytaine
%              linear damping inputs; retain commented per-vessel overrides.
%   2026-10-07 Added generic editable defaults for custom vessel MAT-files.
%   2026-10-09 Use initial JONSWAP parameters within the recommended range.

if nargin < 2
    useGenericDefaults = false;
end
if ~isscalar(useGenericDefaults) || ...
        ~(islogical(useGenericDefaults) || ...
        (isnumeric(useGenericDefaults) && ...
        (useGenericDefaults == 0 || useGenericDefaults == 1)))
    error('useGenericDefaults must be a logical scalar.');
end

% Vessel data
cfg.vessel.matFile = matFile;

% RAO update time (10 Hz)
cfg.RAO_update_period = 0.1; 

% Ocean environment and wave discretization
cfg.environment.Vc = 0.5;
cfg.environment.betaVc = deg2rad(30);
cfg.environment.Hs = 2.0;
cfg.environment.w0 = 1.0;
cfg.environment.beta_wave = deg2rad(140);
cfg.environment.spectrumType = 'JONSWAP';
cfg.environment.spreadingFlag = 1;
cfg.environment.numFreqIntervals = 100;
cfg.environment.numDirections = 24;

% Generic simulation and initial-state defaults. Named vessels override
% these values below, while a custom vessel presents them directly in the GUI.
cfg.simulation.T_final = 180;
cfg.simulation.h = 0.02;
cfg.initial.nu = zeros(6,1);
cfg.initial.eta = zeros(6,1);

% Controller selection and setpoint change. Heading autopilot is the default;
% the supply vessel and semisubmersible override this with DP control.
cfg.control.mode = 'headingAutopilot';
cfg.control.setpointChangeTime = 50;

% Dynamic-positioning controller
cfg.control.dp.eta_ref = [0; 0; 0];
cfg.control.dp.eta_ref_after = [0; 0; deg2rad(40)];
cfg.control.dp.wn = [0.1 0.1 0.3];
cfg.control.dp.zeta = [1 1 1];
cfg.control.dp.T_f = 30;

% Heading autopilot
cfg.control.heading.psi_ref = 0;
cfg.control.heading.psi_ref_after = deg2rad(30);
cfg.control.heading.r_max = deg2rad(2.0);
cfg.control.heading.wn = 1.0;
cfg.control.heading.zeta = 1.0;
cfg.control.heading.tauX = 0;

% Do not assume nonlinear damping for a newly developed vessel.
cfg.damping.nonlinear_456 = [0 0 0];

% Vessel-specific overrides
if ~useGenericDefaults
    switch lower(matFile)
    case 'supply.mat'
        cfg.simulation.T_final = 600;
        cfg.simulation.h = 0.05;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = [10, 10, 0, deg2rad(5), 0, 0]';
        cfg.control.mode = 'DPsystem';
        cfg.control.dp.wn = [0.3 0.3 0.9];
        cfg.damping.kappa_126 = [0.05 0.05 0.05]; % Bv is 5% of B_eq-values DOFs 1,2,6
        cfg.damping.delta_zeta_345 = [0 0.1 0];   % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 0 0];      % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 180e3;

    case 's175.mat'
        cfg.simulation.T_final = 600;
        cfg.simulation.h = 0.05;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = [10, 10, 0, deg2rad(5), 0, 0]';
        cfg.control.dp.wn = [0.3 0.3 0.9];
        cfg.damping.kappa_126 = [0.05 0.05 0.05]; % Bv is 5% of B_eq-values DOFs 1,2,6
        cfg.damping.delta_zeta_345 = [0 0.2 0];   % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 0 0];      % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 120e3;

    case 'tanker.mat'
        cfg.simulation.T_final = 600;
        cfg.simulation.h = 0.05;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = zeros(6,1);
        cfg.damping.kappa_126 = [0.05 0.05 0.05]; % Bv is 5% of B_eq-values DOFs 1,2,6
        cfg.damping.delta_zeta_345 = [0 0.1 0];   % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 0 0];      % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 120e3;

    case 'semisub.mat'
        cfg.simulation.T_final = 500;
        cfg.simulation.h = 0.05;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = zeros(6,1);
        cfg.control.mode = 'DPsystem';
        cfg.damping.kappa_126 = [0.05 0.05 0.05];  % Bv is 5% of B_eq-values DOFs 1,2,6
        cfg.damping.delta_zeta_345 = [0 0.1 0];    % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 0 0];       % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 120e3;

    case 'fpso.mat'
        cfg.simulation.T_final = 600;
        cfg.simulation.h = 0.05;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = zeros(6,1);
        cfg.damping.kappa_126 = [0.05 0.05 0.05]; % Bv is 5% of B_eq-values DOFs 1,2,6
        cfg.damping.delta_zeta_345 = [0 0.1 0];   % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 0 0];      % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 120e3;

    case 'testship.mat'
        cfg.simulation.T_final = 180;
        cfg.simulation.h = 0.02;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = [10, 10, 0, deg2rad(5), 0, 0]';
        cfg.control.dp.wn = [0.3 0.3 0.9];
        % Optional overrides of the MSS-Capytaine values exported from config.json:
        % cfg.damping.kappa_126 = [0.05 0.05 0.05]; % Bv is 5% of B_eq-values DOFs 1,2,6
        % cfg.damping.delta_zeta_345 = [0 0.2 0];   % Increase damping ratios DOFs 3,4,5
        cfg.damping.nonlinear_456 = [5 5 0];        % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.wn = 1.0;               % M = Mq * (1 + kappa5 * |q|)
        cfg.control.heading.zeta = 1.0;
        cfg.control.heading.tauX = 10e3;

    case 'lauv_marie.mat'
        cfg.simulation.T_final = 180;
        cfg.simulation.h = 0.02;
        cfg.environment.Hs = 0.5;
        cfg.environment.w0 = 2.0;
        cfg.initial.nu = zeros(6,1);
        cfg.initial.eta = [0, 0, 5, 0, 0, 0]';
        cfg.control.dp.wn = [0.5 0.5 1.5];
        % Optional overrides of the MSS-Capytaine values exported from config.json:
        % cfg.damping.T_1236 = [50 5 5 5];       % Time constants in DOFs 1,2,3,6
        % cfg.damping.delta_zeta_45 = [0.2 0.2]; % Increase damping ratios DOFs 4,5
        cfg.damping.nonlinear_456 = [5 5 0];     % K = Kp * (1 + kappa4 * |p|)
        cfg.control.heading.tauX = 20;           % M = Mq * (1 + kappa5 * |q|)
        cfg.control.heading.wn = 1.5;
        cfg.control.heading.zeta = 1.0;

        otherwise
            error('Unsupported hydrodynamic model: %s', matFile);
    end
end

cfg = useStoredPowerBasedDamping(cfg);

if useGenericDefaults
    cfg = useStoredSubmergenceDepth(cfg);
end

end

function cfg = useStoredSubmergenceDepth(cfg)
% Initialize a custom submerged vehicle at its exported operating depth.

matPath = resolveHydroVesselFile(cfg.vessel.matFile);
data = load(matPath, 'vessel');
if ~isfield(data, 'vessel') || ~isfield(data.vessel, 'main')
    return
end

main = data.vessel.main;
if isfield(main, 'submerged') && isscalar(main.submerged) && ...
        logical(main.submerged) && isfield(main, 'submergenceDepth') && ...
        isscalar(main.submergenceDepth) && isfinite(main.submergenceDepth)
    cfg.initial.eta(3) = main.submergenceDepth;
end

end

function cfg = useStoredPowerBasedDamping(cfg)
% Use explicit cfg values as overrides; otherwise load the MAT-file defaults.

hasKappa = isfield(cfg.damping, 'kappa_126');
hasDelta345 = isfield(cfg.damping, 'delta_zeta_345');
hasTimeConstants = isfield(cfg.damping, 'T_1236');
hasDelta45 = isfield(cfg.damping, 'delta_zeta_45');

if hasKappa || hasDelta345
    if ~hasKappa || ~hasDelta345 || hasTimeConstants || hasDelta45
        error(['A floating-vessel damping override requires both ' ...
            'kappa_126 and delta_zeta_345, and no submerged inputs.']);
    end
    return
end

if hasTimeConstants || hasDelta45
    if ~hasTimeConstants || ~hasDelta45
        error(['A submerged-vehicle damping override requires both ' ...
            'T_1236 and delta_zeta_45.']);
    end
    return
end

matPath = resolveHydroVesselFile(cfg.vessel.matFile);
data = load(matPath, 'vessel');
if ~isfield(data, 'vessel') || ...
        ~isfield(data.vessel, 'powerBased')
    error(['%s does not contain power-based damping inputs. Add an ' ...
        'explicit damping override in hydroVesselConfig.m.'], ...
        cfg.vessel.matFile);
end

powerBased = data.vessel.powerBased;
if isfield(powerBased, 'T_1236') && ...
        isfield(powerBased, 'delta_zeta_45')
    cfg.damping.T_1236 = powerBased.T_1236;
    cfg.damping.delta_zeta_45 = powerBased.delta_zeta_45;
elseif isfield(powerBased, 'kappa_126') && ...
        isfield(powerBased, 'delta_zeta_345')
    cfg.damping.kappa_126 = powerBased.kappa_126;
    cfg.damping.delta_zeta_345 = powerBased.delta_zeta_345;
else
    error(['%s has incomplete power-based damping inputs. Add an ' ...
        'explicit damping override in hydroVesselConfig.m.'], ...
        cfg.vessel.matFile);
end

end
