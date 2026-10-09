function [vessel, cfg] = guiSIMhydroVessel()
% [vessel,cfg] = guiSIMhydroVessel() displays the hydrodynamic vessel
% simulation options. Selecting a vessel loads its defaults from
% hydroVesselConfig. Pressing OK returns the vessel data and the edited
% configuration structure. Closing the window returns empty outputs.
%
% Author:    Thor I. Fossen
% Date:      2026-10-04
% Revisions:
%   2026-10-04 Added user-facing controller names while preserving the
%              internal DPsystem and headingAutopilot mode flags.
%   2026-10-07 Added loading of custom MSS-compatible vessel MAT-files.

defaultMatFile = 'testShip.mat';
cfg = hydroVesselConfig(defaultMatFile);
vessel = [];
accepted = false;
selectedMatFile = defaultMatFile;
useGenericDefaults = false;
previousSelectedRadio = [];
customStartPath = fileparts(resolveHydroVesselFile(defaultMatFile));

f = figure('Position', [120, 60, 1200, 820], ...
    'Name', 'Hydrodynamic Vessel Simulation Options', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'Resize', 'off', ...
    'WindowStyle', 'modal', ...
    'CloseRequestFcn', @onCancel);

% Vessel selection
bgVessel = uibuttongroup('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [15 560 365 245], ...
    'Title', 'Marine craft', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

addText(bgVessel, 'ShipX (commercial)', [10 185 160 24], 'bold');
addRadio(bgVessel, 'Supply vessel', 'supply.mat', [20 155 150 26], 0);
addRadio(bgVessel, 'S175 container ship', 's175.mat', [20 127 160 26], 0);

addText(bgVessel, 'WAMIT (commercial)', [185 185 165 24], 'bold');
addRadio(bgVessel, 'Tanker', 'tanker.mat', [195 155 145 26], 0);
addRadio(bgVessel, 'FPSO', 'fpso.mat', [195 127 145 26], 0);
addRadio(bgVessel, 'Semisubmersible', 'semisub.mat', [195 99 150 26], 0);

addText(bgVessel, 'Capytaine (open source)', [10 65 190 24], 'bold');
radioTestShip = addRadio(bgVessel, 'Test ship', defaultMatFile, ...
    [20 35 95 26], 1);
addRadio(bgVessel, 'Light AUV (LAUV) Marie', 'LAUV_marie.mat', ...
    [120 35 220 26], 0);
radioCustom = addRadio(bgVessel, 'Custom vessel', '__custom__', ...
    [20 5 155 26], 0);
uicontrol('Parent', bgVessel, ...
    'Style', 'pushbutton', ...
    'String', '1. Load MAT file...', ...
    'Position', [190 4 150 28], ...
    'FontSize', 10, ...
    'Callback', @onBrowseCustom);
previousSelectedRadio = radioTestShip;

% Environmental parameters
pEnvironment = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [395 455 360 350], ...
    'Title', 'Environment', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

addText(pEnvironment, 'Wave spectrum', [15 295 190 22], 'normal');
popupSpectrum = uicontrol('Parent', pEnvironment, ...
    'Style', 'popupmenu', ...
    'String', {'Modified PM', 'JONSWAP', 'Torsethaugen'}, ...
    'Position', [220 295 120 25], ...
    'FontSize', 10);

checkboxSpreading = uicontrol('Parent', pEnvironment, ...
    'Style', 'checkbox', ...
    'String', 'Enable directional spreading', ...
    'Position', [15 264 240 25], ...
    'FontSize', 10);

editHs = addEdit(pEnvironment, 'Significant wave height Hs (m)', 230, 220, 110);
editW0 = addEdit(pEnvironment, 'Peak frequency w0 (rad/s)', 198, 220, 110);
editBetaWave = addEdit(pEnvironment, 'Wave direction (deg)', 166, 220, 110);
editVc = addEdit(pEnvironment, 'Current speed Vc (m/s)', 134, 220, 110);
editBetaVc = addEdit(pEnvironment, 'Current direction (deg)', 102, 220, 110);
editNumFreq = addEdit(pEnvironment, 'Frequency intervals (> 50)', 70, 220, 110);
editNumDirections = addEdit(pEnvironment, 'Wave directions (> 15)', 38, 220, 110);

% Simulation and initial state
pSimulation = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [770 455 415 350], ...
    'Title', 'Simulation and initial state', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

editTFinal = addEdit(pSimulation, 'Simulation time (s)', 295, 195, 195);
editH = addEdit(pSimulation, 'Integration step h (s)', 263, 195, 195);
addText(pSimulation, ...
    sprintf('RAO update period: %g s (fixed)', cfg.RAO_update_period), ...
    [15 231 385 22], 'normal');
editChangeTime = addEdit(pSimulation, 'Setpoint change time (s)', 199, 195, 195);
editEta0 = addEdit(pSimulation, ...
    'eta = [xn yn zn theta psi]', 150, 195, 195);
addText(pSimulation, 'Positions in m; angles in deg', ...
    [15 125 195 20], 'normal');
editNu0 = addEdit(pSimulation, ...
    'nu = [u v w p q r]', 92, 195, 195);
addText(pSimulation, 'Linear velocities in m/s; angular rates in rad/s', ...
    [15 60 385 22], 'normal');

% Controller parameters
pControl = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [15 5 740 440], ...
    'Title', 'Control system', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

addText(pControl, 'Control system', [15 365 100 22], 'bold');
addText(pControl, ...
    'The displayed controller gains are editable starting values.', ...
    [315 365 405 22], 'normal');
controlLabels = {'DP control system', 'Heading autopilot'};
controlModes = {'DPsystem', 'headingAutopilot'};
popupControl = uicontrol('Parent', pControl, ...
    'Style', 'popupmenu', ...
    'String', controlLabels, ...
    'Position', [120 365 180 25], ...
    'FontSize', 10);

pDP = uipanel('Parent', pControl, ...
    'Units', 'pixels', ...
    'Position', [10 10 350 340], ...
    'Title', 'Dynamic positioning', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

editEtaRef = addEdit(pDP, 'Initial eta_ref [x y psi]', 285, 175, 150);
editEtaRefAfter = addEdit(pDP, 'After change [x y psi]', 245, 175, 150);
addText(pDP, 'Positions in m; heading in deg', [15 218 315 20], 'normal');
editDPWn = addEdit(pDP, 'wn [x y psi] (rad/s)', 180, 175, 150);
editDPZeta = addEdit(pDP, 'zeta [x y psi]', 140, 175, 150);
editTf = addEdit(pDP, 'Setpoint filter T_f (s)', 100, 175, 150);

pHeading = uipanel('Parent', pControl, ...
    'Units', 'pixels', ...
    'Position', [375 10 350 340], ...
    'Title', 'Heading autopilot', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

editPsiRef = addEdit(pHeading, 'Initial heading (deg)', 285, 175, 150);
editPsiRefAfter = addEdit(pHeading, 'Heading after change (deg)', 245, 175, 150);
editRMax = addEdit(pHeading, 'Maximum yaw rate (deg/s)', 205, 175, 150);
editHeadingWn = addEdit(pHeading, 'Natural frequency (rad/s)', 165, 175, 150);
editHeadingZeta = addEdit(pHeading, 'Relative damping ratio', 125, 175, 150);
editTauX = addEdit(pHeading, 'Constant surge force (N)', 85, 175, 150);

% Damping parameters
pDamping = uipanel('Parent', f, ...
    'Units', 'pixels', ...
    'Position', [770 175 415 270], ...
    'Title', 'Viscous damping', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

editAperiodicDamping = addEdit(pDamping, 'Aperiodic damping parameters', ...
    205, 220, 170);
textAperiodicDamping = addText(pDamping, 'kappa_126: [surge sway yaw]', ...
    [15 180 250 20], 'normal');
editDeltaZeta = addEdit(pDamping, 'Damping ratio increase', ...
    140, 220, 170);
textDeltaZeta = addText(pDamping, 'delta_zeta_345: [heave roll pitch]', ...
    [15 115 250 20], 'normal');
editNonlinear456 = addEdit(pDamping, 'Quadratic damping factors', ...
    75, 220, 170);
addText(pDamping, 'kappa_456: [roll pitch yaw]', ...
    [15 50 250 20], 'normal');
addText(pDamping, ...
    'Vessel defaults are reloaded when another craft is selected.', ...
    [15 10 380 25], 'normal');

% Current selection and next action
textStatus = addText(f, ...
    {'1. Selected: Test ship    2. Choose gains and settings', ...
    '3. Click Run simulation'}, ...
    [770 128 415 38], 'bold');

% Dialog buttons
uicontrol('Parent', f, ...
    'Style', 'pushbutton', ...
    'String', '3. Run simulation', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'Position', [865 82 145 42], ...
    'Callback', @onOK);

uicontrol('Parent', f, ...
    'Style', 'pushbutton', ...
    'String', 'Cancel', ...
    'FontSize', 11, ...
    'Position', [1025 85 100 38], ...
    'Callback', @onCancel);

populateControls(cfg);
uiwait(f);

if ishandle(f)
    delete(f);
end

if ~accepted
    vessel = [];
    cfg = [];
end

    function radio = addRadio(parent, label, tag, position, value)
        radio = uicontrol('Parent', parent, ...
            'Style', 'radiobutton', ...
            'String', label, ...
            'Tag', tag, ...
            'Value', value, ...
            'Position', position, ...
            'FontSize', 10, ...
            'Callback', @onVesselChanged);
    end

    function onVesselChanged(source, ~)
        if get(source, 'Value') ~= 1
            return
        end

        if strcmp(get(source, 'Tag'), '__custom__')
            chooseCustomFile(source);
            return
        end

        try
            newMatFile = get(source, 'Tag');
            newCfg = hydroVesselConfig(newMatFile);
            populateControls(newCfg);
            selectedMatFile = newMatFile;
            useGenericDefaults = false;
            previousSelectedRadio = source;
            cfg = newCfg;
            set(f, 'Name', 'Hydrodynamic Vessel Simulation Options');
            set(textStatus, 'String', ...
                {['1. Selected: ', get(source, 'String'), ...
                '    2. Choose gains and settings'], ...
                '3. Click Run simulation'});
        catch exception
            restorePreviousSelection(source);
            errordlg(exception.message, 'Unable to load vessel', 'modal');
        end
    end

    function onBrowseCustom(~, ~)
        chooseCustomFile(radioCustom);
    end

    function chooseCustomFile(source)
        [fileName, pathName] = uigetfile( ...
            {'*.mat', 'MAT-files (*.mat)'}, ...
            'Load custom hydrodynamic vessel', ...
            fullfile(customStartPath, '*.mat'));
        if isequal(fileName, 0)
            restorePreviousSelection(source);
            return
        end

        try
            matPath = fullfile(pathName, fileName);
            loadHydroVesselFile(matPath);
            newCfg = hydroVesselConfig(matPath, true);
            populateControls(newCfg);

            selectRadio(radioCustom);
            set(radioCustom, 'String', ['Custom: ', fileName]);
            selectedMatFile = matPath;
            useGenericDefaults = true;
            previousSelectedRadio = radioCustom;
            customStartPath = pathName;
            cfg = newCfg;
            set(f, 'Name', ...
                ['Hydrodynamic Vessel Simulation Options - ', fileName]);
            set(textStatus, 'String', ...
                {['1. Loaded: ', fileName, ...
                '    2. Choose gains and settings'], ...
                '3. Click Run simulation'});
        catch exception
            restorePreviousSelection(source);
            errordlg(exception.message, 'Unable to load vessel', 'modal');
        end
    end

    function restorePreviousSelection(source)
        if ~isempty(previousSelectedRadio) && ...
                ishandle(previousSelectedRadio)
            selectRadio(previousSelectedRadio);
        elseif ishandle(source)
            set(source, 'Value', 0);
        end
    end

    function selectRadio(target)
        radios = findobj(bgVessel, 'Style', 'radiobutton');
        set(radios, 'Value', 0);
        set(target, 'Value', 1);
    end

    function populateControls(defaults)
        spectrumNames = get(popupSpectrum, 'String');
        spectrumIndex = find(strcmp(spectrumNames, ...
            defaults.environment.spectrumType), 1);
        if isempty(spectrumIndex)
            spectrumIndex = 1;
        end
        set(popupSpectrum, 'Value', spectrumIndex);
        set(checkboxSpreading, 'Value', defaults.environment.spreadingFlag);
        set(editHs, 'String', scalarText(defaults.environment.Hs));
        set(editW0, 'String', scalarText(defaults.environment.w0));
        set(editBetaWave, 'String', ...
            scalarText(rad2deg(defaults.environment.beta_wave)));
        set(editVc, 'String', scalarText(defaults.environment.Vc));
        set(editBetaVc, 'String', ...
            scalarText(rad2deg(defaults.environment.betaVc)));
        set(editNumFreq, 'String', ...
            scalarText(defaults.environment.numFreqIntervals));
        set(editNumDirections, 'String', ...
            scalarText(defaults.environment.numDirections));

        set(editTFinal, 'String', scalarText(defaults.simulation.T_final));
        set(editH, 'String', scalarText(defaults.simulation.h));
        set(editChangeTime, 'String', ...
            scalarText(defaults.control.setpointChangeTime));

        eta0 = defaults.initial.eta(:)';
        eta0(4:6) = rad2deg(eta0(4:6));
        set(editEta0, 'String', vectorText(eta0));
        set(editNu0, 'String', vectorText(defaults.initial.nu));

        controlIndex = find(strcmp(controlModes, defaults.control.mode), 1);
        if isempty(controlIndex)
            controlIndex = 1;
        end
        set(popupControl, 'Value', controlIndex);

        etaRef = defaults.control.dp.eta_ref(:)';
        etaRef(3) = rad2deg(etaRef(3));
        etaRefAfter = defaults.control.dp.eta_ref_after(:)';
        etaRefAfter(3) = rad2deg(etaRefAfter(3));
        set(editEtaRef, 'String', vectorText(etaRef));
        set(editEtaRefAfter, 'String', vectorText(etaRefAfter));
        set(editDPWn, 'String', vectorText(defaults.control.dp.wn));
        set(editDPZeta, 'String', vectorText(defaults.control.dp.zeta));
        set(editTf, 'String', scalarText(defaults.control.dp.T_f));

        set(editPsiRef, 'String', ...
            scalarText(rad2deg(defaults.control.heading.psi_ref)));
        set(editPsiRefAfter, 'String', ...
            scalarText(rad2deg(defaults.control.heading.psi_ref_after)));
        set(editRMax, 'String', ...
            scalarText(rad2deg(defaults.control.heading.r_max)));
        set(editHeadingWn, 'String', ...
            scalarText(defaults.control.heading.wn));
        set(editHeadingZeta, 'String', ...
            scalarText(defaults.control.heading.zeta));
        set(editTauX, 'String', ...
            scalarText(defaults.control.heading.tauX));

        if isfield(defaults.damping, 'T_1236')
            set(editAperiodicDamping, 'String', ...
                vectorText(defaults.damping.T_1236));
            set(textAperiodicDamping, 'String', ...
                'T_1236 [s]: [surge sway heave yaw]');
            set(editDeltaZeta, 'String', ...
                vectorText(defaults.damping.delta_zeta_45));
            set(textDeltaZeta, 'String', ...
                'delta_zeta_45: [roll pitch]');
        else
            set(editAperiodicDamping, 'String', ...
                vectorText(defaults.damping.kappa_126));
            set(textAperiodicDamping, 'String', ...
                'kappa_126: [surge sway yaw]');
            set(editDeltaZeta, 'String', ...
                vectorText(defaults.damping.delta_zeta_345));
            set(textDeltaZeta, 'String', ...
                'delta_zeta_345: [heave roll pitch]');
        end
        set(editNonlinear456, 'String', ...
            vectorText(defaults.damping.nonlinear_456));
    end

    function onOK(~, ~)
        try
            matFile = selectedMatFile;
            newCfg = hydroVesselConfig(matFile, useGenericDefaults);

            spectrumNames = get(popupSpectrum, 'String');
            newCfg.environment.spectrumType = ...
                spectrumNames{get(popupSpectrum, 'Value')};
            newCfg.environment.spreadingFlag = ...
                get(checkboxSpreading, 'Value');
            newCfg.environment.Hs = readScalar(editHs, 'Hs');
            newCfg.environment.w0 = readScalar(editW0, 'w0');
            newCfg.environment.beta_wave = deg2rad( ...
                readScalar(editBetaWave, 'wave direction'));
            newCfg.environment.Vc = readScalar(editVc, 'Vc');
            newCfg.environment.betaVc = deg2rad( ...
                readScalar(editBetaVc, 'current direction'));
            newCfg.environment.numFreqIntervals = ...
                readScalar(editNumFreq, 'frequency intervals');
            newCfg.environment.numDirections = ...
                readScalar(editNumDirections, 'wave directions');

            newCfg.simulation.T_final = ...
                readScalar(editTFinal, 'simulation time');
            newCfg.simulation.h = readScalar(editH, 'integration step');
            newCfg.control.setpointChangeTime = ...
                readScalar(editChangeTime, 'setpoint change time');

            eta0 = readVector(editEta0, 'initial eta', 6);
            eta0(4:6) = deg2rad(eta0(4:6));
            newCfg.initial.eta = eta0(:);
            newCfg.initial.nu = readVector(editNu0, 'initial nu', 6)';

            newCfg.control.mode = ...
                controlModes{get(popupControl, 'Value')};

            etaRef = readVector(editEtaRef, 'initial DP reference', 3);
            etaRef(3) = deg2rad(etaRef(3));
            newCfg.control.dp.eta_ref = etaRef(:);
            etaRefAfter = readVector(editEtaRefAfter, ...
                'DP reference after change', 3);
            etaRefAfter(3) = deg2rad(etaRefAfter(3));
            newCfg.control.dp.eta_ref_after = etaRefAfter(:);
            newCfg.control.dp.wn = readVector(editDPWn, 'DP wn', 3);
            newCfg.control.dp.zeta = ...
                readVector(editDPZeta, 'DP zeta', 3);
            newCfg.control.dp.T_f = readScalar(editTf, 'DP filter time');

            newCfg.control.heading.psi_ref = deg2rad( ...
                readScalar(editPsiRef, 'initial heading reference'));
            newCfg.control.heading.psi_ref_after = deg2rad( ...
                readScalar(editPsiRefAfter, 'heading reference after change'));
            newCfg.control.heading.r_max = deg2rad( ...
                readScalar(editRMax, 'maximum yaw rate'));
            newCfg.control.heading.wn = ...
                readScalar(editHeadingWn, 'heading natural frequency');
            newCfg.control.heading.zeta = ...
                readScalar(editHeadingZeta, 'heading damping ratio');
            newCfg.control.heading.tauX = ...
                readScalar(editTauX, 'constant surge force');

            if isfield(newCfg.damping, 'T_1236')
                newCfg.damping.T_1236 = ...
                    readVector(editAperiodicDamping, 'T_1236', 4);
                newCfg.damping.delta_zeta_45 = ...
                    readVector(editDeltaZeta, 'delta_zeta_45', 2);
            else
                newCfg.damping.kappa_126 = ...
                    readVector(editAperiodicDamping, 'kappa_126', 3);
                newCfg.damping.delta_zeta_345 = ...
                    readVector(editDeltaZeta, 'delta_zeta_345', 3);
            end
            newCfg.damping.nonlinear_456 = ...
                readVector(editNonlinear456, 'nonlinear kappa_456', 3);

            validateConfiguration(newCfg);

            matPath = resolveHydroVesselFile(matFile);
            data = loadHydroVesselFile(matPath);

            vessel = data.vessel;
            vessel.matFile = matFile;
            vessel.kappa_4 = newCfg.damping.nonlinear_456(1);
            vessel.kappa_5 = newCfg.damping.nonlinear_456(2);
            vessel.kappa_6 = newCfg.damping.nonlinear_456(3);
            cfg = newCfg;
            accepted = true;
            uiresume(f);

        catch exception
            errordlg(exception.message, 'Invalid configuration', 'modal');
        end
    end

    function onCancel(~, ~)
        accepted = false;
        uiresume(f);
    end

end

function data = loadHydroVesselFile(matPath)
% Load and minimally validate a hydrodynamic vessel MAT-file.

data = load(matPath, 'vessel');
if ~isfield(data, 'vessel') || ~isstruct(data.vessel) || ...
        ~isscalar(data.vessel)
    error('%s does not contain a vessel structure.', matPath);
end

% powerBased is optional here: legacy ShipX and WAMIT files build it when
% computeManeuveringModel runs, while configurations that need stored damping
% inputs validate it in hydroVesselConfig.
requiredFields = {'main', 'MRB', 'A', 'B', 'C', 'forceRAO', ...
    'freqs', 'headings', 'velocities'};
missingFields = requiredFields(~isfield(data.vessel, requiredFields));
if ~isempty(missingFields)
    error('%s is missing vessel field(s): %s.', matPath, ...
        strjoin(missingFields, ', '));
end

end

function handle = addEdit(parent, label, yPosition, labelWidth, editWidth)
addText(parent, label, [15 yPosition labelWidth 22], 'normal');
handle = uicontrol('Parent', parent, ...
    'Style', 'edit', ...
    'String', '', ...
    'HorizontalAlignment', 'left', ...
    'Position', [labelWidth + 15 yPosition editWidth 24], ...
    'FontSize', 10);
end

function handle = addText(parent, label, position, fontWeight)
handle = uicontrol('Parent', parent, ...
    'Style', 'text', ...
    'String', label, ...
    'HorizontalAlignment', 'left', ...
    'Position', position, ...
    'FontSize', 10, ...
    'FontWeight', fontWeight);
end

function value = readScalar(handle, label)
value = str2double(get(handle, 'String'));
if ~isscalar(value) || ~isfinite(value)
    error('%s must be a finite scalar.', label);
end
end

function value = readVector(handle, label, expectedLength)
inputText = get(handle, 'String');
inputText = regexprep(inputText, '[\[\],;]', ' ');
value = sscanf(inputText, '%f')';
if numel(value) ~= expectedLength || any(~isfinite(value))
    error('%s must contain %d finite values.', label, expectedLength);
end
end

function validateConfiguration(cfg)
if cfg.simulation.h <= 0 || cfg.simulation.h > 0.1
    error('The integration step h must be positive and no greater than 0.1 s.');
end
if cfg.simulation.T_final <= 0
    error('The simulation time must be positive.');
end
if cfg.control.setpointChangeTime < 0
    error('The setpoint change time cannot be negative.');
end
if cfg.environment.Hs <= 0 || cfg.environment.w0 <= 0
    error('Hs and w0 must be positive.');
end
if cfg.environment.Vc < 0
    error('The current speed cannot be negative.');
end
if cfg.environment.numFreqIntervals <= 50 || ...
        cfg.environment.numFreqIntervals ~= ...
        round(cfg.environment.numFreqIntervals)
    error('The number of frequency intervals must be an integer above 50.');
end
if cfg.environment.numDirections <= 15 || ...
        cfg.environment.numDirections ~= round(cfg.environment.numDirections)
    error('The number of wave directions must be an integer above 15.');
end
if any(cfg.control.dp.wn <= 0) || any(cfg.control.dp.zeta < 0) || ...
        cfg.control.dp.T_f <= 0
    error('The DP tuning parameters must be nonnegative, with wn and T_f positive.');
end
if cfg.control.heading.r_max <= 0 || cfg.control.heading.wn <= 0 || ...
        cfg.control.heading.zeta < 0
    error('The heading-controller parameters are invalid.');
end
if isfield(cfg.damping, 'T_1236')
    invalidLinearDamping = any(cfg.damping.T_1236 <= 0) || ...
        any(cfg.damping.delta_zeta_45 < 0);
else
    invalidLinearDamping = any(cfg.damping.kappa_126 < 0) || ...
        any(cfg.damping.delta_zeta_345 < 0);
end
if invalidLinearDamping || any(cfg.damping.nonlinear_456 < 0)
    error(['Time constants must be positive; damping increments and ' ...
        'nonlinear damping factors must be nonnegative.']);
end
end

function textValue = scalarText(value)
textValue = sprintf('%.12g', value);
end

function textValue = vectorText(value)
textValue = mat2str(value(:)', 12);
end
