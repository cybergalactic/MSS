function SIMhydroVessel()
% The function simulates a selected ShipX, WAMIT, or Capytaine vessel model using
% power-based constant hydrodynamic matrices and 1st-order force RAOs. The 
% control system can be chosen as a nonlinear MIMO PID controller for DP or a
% heading autopilot for transit. Vessel-specific defaults are defined by
% hydroVesselConfig.m and can be edited in guiSIMhydroVessel.m before simulation.
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-10-04
% Revisions:

clear PIDnonlinearMIMO headingAutopilot; % Clear persistent state variables 
clearvars;                               % Clear all other variables
close all;                               % Close all windows

% Load the selected vessel and its config. parameters; see hydroVesselConfig.m
[vessel, cfg] = guiSIMhydroVessel(); 
disp(['Marine craft: ', vessel.matFile]); % Display the vessel name    

% ------------------------------------------------------------------------------
% Wave spectrum and force RAO initializations
% ------------------------------------------------------------------------------
environment = cfg.environment;
rng(1); % Reproducible stochastic wave realization
[environment, waveModel] = waveInitialization(vessel, environment);

% ------------------------------------------------------------------------------
% Compute power-based equivalent matrices and viscous damping:
%   Fossen, T. I. (2025). Maneuvering Coefficient Estimation from Frequency-
%   Dependent Added Mass and Damping: A Power-Based Approach. Ocean Engineering, 
%   341, 122494.
% ------------------------------------------------------------------------------
if isfield(cfg.damping, 'T_1236')
    % Use time constants in surge, sway, heave, and yaw for underwater vehicles
    aperiodicDamping = cfg.damping.T_1236;
    delta_zeta = cfg.damping.delta_zeta_45;
else
    % Use percentage damping increase in surge, sway, and yaw for surface craft
    aperiodicDamping = cfg.damping.kappa_126;
    delta_zeta = cfg.damping.delta_zeta_345;
end
vessel = computeManeuveringModel(vessel, environment.w0, ...
    aperiodicDamping, delta_zeta, 0);

% ------------------------------------------------------------------------------
% Time vector initialization
% ------------------------------------------------------------------------------
t = 0:cfg.simulation.h:cfg.simulation.T_final; % Time vector, sampling time h
nextRAOtime = 0;                               % Next RAO update time
nTimeSteps = length(t);                        % Number of time steps

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
simdata = zeros(nTimeSteps, 25); % Pre-allocate matrix for efficiency
x = [cfg.initial.nu; cfg.initial.eta]; % Initial state vector

for i = 1:nTimeSteps

    % Measurements
    nu = x(1:6);
    eta = x(7:12);

    % Control logic based on the elapsed simulation time
    if t(i) > cfg.control.setpointChangeTime
        cfg.control.dp.eta_ref = cfg.control.dp.eta_ref_after; % DP
        cfg.control.heading.psi_ref = cfg.control.heading.psi_ref_after; % Autopilot
    end

    switch cfg.control.mode

        case 'DPsystem'

            % Nonlinear MIMO PID controller for dynamic positioning (DP)
            tau = PIDnonlinearMIMO(...
                eta, nu, cfg.control.dp.eta_ref, diag(diag(vessel.M)), ...
                diag(cfg.control.dp.wn), diag(cfg.control.dp.zeta), ...
                cfg.control.dp.T_f, cfg.simulation.h);

        case 'headingAutopilot'

            % Heading autopilot PID controller
            tauN = headingAutopilot(eta(6), nu(6), ...
                cfg.control.heading.psi_ref, cfg.control.heading.r_max, ...
                vessel.M(6,6), cfg.control.heading.wn, ...
                cfg.control.heading.zeta, cfg.simulation.h);
            tau = [cfg.control.heading.tauX; 0; 0; 0; 0; tauN];

        otherwise
            error('Unsupported control system: %s', cfg.control.mode);

    end

    % 6-DOF generalized wave forces (compute RAO only every 0.1 second)
    if t(i) >= nextRAOtime
        U = sqrt(nu(1)^2 + nu(2)^2);
        [tau_wave1, waveElevation] = waveForceRAO(t(i), waveModel, ...
            U, eta(6), environment.beta_wave);

        nextRAOtime = nextRAOtime + cfg.RAO_update_period;
    end

    % Save data for post-processing and plotting
    simdata(i, :) = [eta; nu; tau; tau_wave1; waveElevation];

    % RK4 integration (k+1)
    tau_ext = tau + tau_wave1; % Sum of external forces
    x = rk4(@hydroVessel, cfg.simulation.h, x, tau_ext, vessel, ...
        environment.Vc, environment.betaVc);

end

% ------------------------------------------------------------------------------
% Plot the simulation data and display main characteristics
% ------------------------------------------------------------------------------
plotSIMhydroVessel(t, simdata, environment); 
modalFromMDG(vessel.M, vessel.D, vessel.G, 1);
displayHydroData(vessel); 

end
