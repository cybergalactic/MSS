% This script simulates the response of a vessel in waves using Cummins (1962) 
% equation and an equivalent maneuvering model (Fossen, 2025). It calculates 
% the wave-induced forces using 1st-order force RAOs, solves the full hydrodynamic 
% model including memory effects, and compares it with a simplified approximation 
% using constant equivalent added mass and damping matrices (A_eq, B_eq) according 
% according to: 
%
%   A_eq = ∫ A_U(ω) S_N(ω) dω
%   B_eq = ∫ B_U(ω) S_N(ω) dω 
%
% where
%
%   S(ω)                 : Wave energy spectrum
%   S_N(ω) = S(ω) / m_0  : Normalized wave spectrum
%   m_0 = ∫ S(ω) dω      : Zero spectral moment
%   A_U(ω) and B_U(ω)    : Added mass and damping matrices at speed U
%  
% The integrals are evaluated numerically using trapezoidal integration.
% This ensures that the kinetic energy and power dissipation properties 
% of the frequency-dependent system are preserved in the equivalent 
% constant matrices.
%
% References:
%   Cummins, W. E. (1962). The impulse response function and ship motions. 
%   Schiffstechnik, 9(47), 101–109.
%
%   Fossen, T. I. (2025). Maneuvering Coefficient Estimation from Frequency-
%   Dependent Added Mass and Damping: A Power-Based Approach. 
%   Ocean Engineering, 341, 122494.
%
% Author:    Thor I. Fossen
% Date:      2025-03-10

clear waveForceRAO; % Clear persistent RAO tables
clearvars; 
close all;
rng(1); % Set random generator seed to 1 when generating stochastic waves

%% SIMULATOR CONFIGURATION
h  = 0.05; % Sampling time (s)
T_final = 120; % Final simulation time (s)
vesselChoice = 1; % Choose vessel type 1, 2, 3

switch vesselChoice
    case 1
        load supply; 
        vesselType = 'Supply Vessel';
        U = vessel.velocities(1); % Zero speed (m/s)
    case 2
        load s175; 
        vesselType = 'S175 Container Ship';
        U = vessel.velocities(3); % Non-zero speed (m/s)
    case 3
        load tanker; 
        vesselType = 'Tanker';
        U = vessel.velocities(1); % Zero speed (m/s)        
end
fprintf('Loaded the %s at %.2f m/s\n', vesselType, U);

psi = 0; % Heading angle (rad)
beta_wave = deg2rad(135); % Wave direction, 0 for following sea, 180 for head sea
maxFreq = 3.0; % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 60; % Number of wave frequency intervals (>50)

% Sea state
Hs = 5; % Significant wave height (m)
omega_p = 0.7;  % Wave spectrum peak frequencies (rad/s)

% First-order force RAOs: Calculate the wave spectrum S(Omega) for each frequency
spectrumType = 'JONSWAP'; 
gamma = 3.3; % Peakedness factor 
Parameter = [Hs, omega_p, gamma]; % Spectrum parameters

% Time vector from 0 to T_final     
t = 0:h:T_final;      
nTimeSteps = length(t);

% Wave spectrum, one direction
omegaMax = vessel.forceRAO.w(end); % Max frequency in RAO dataset

[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumType, ...
    Parameter, numFreqIntervals, omegaMax);

% 6-DOF generalized wave forces using first-order force RAOs
waveData = zeros(nTimeSteps,7); % Pre-allocate table
for i = 1:nTimeSteps
    [tau_wave1, waveElevation] = waveForceRAO(t(i), ...
        S_M, Amp, Omega, mu, vessel, U, psi, beta_wave, numFreqIntervals);
    waveData(i,:) = [tau_wave1' waveElevation];
end

%% Compute Aeq and Beq using simplified normalized wave spectrum
g = 9.81;
omega_p = omega_p - (omega_p^2 / g) * U * cos(beta_wave);
vessel = computeManeuveringModel(vessel, omega_p);

%% Compute Cummins and Maneuvering Model Responses
freqs = vessel.freqs;
nFreqInterp = 200;
freqs_uniform = linspace(min(freqs), max(freqs), nFreqInterp)';

% Initialize storage for all DOFs
eta_cummins = zeros(nTimeSteps,6);  % Displacement (Cummins)
eta_eq = zeros(nTimeSteps,6);       % Displacement (Aeq-Beq)
eta_dot = zeros(nTimeSteps,6);      % Velocity (Cummins)
eta_ddot = zeros(nTimeSteps,6);     % Acceleration (Cummins)
eta_dot_eq = zeros(nTimeSteps,6);   % Velocity (Maneuvering)
A_eq = zeros(6,1);
B_eq = zeros(6,1);
Bv = zeros(6,1);

A_w_all = zeros(length(freqs), 6);
B_w_all = zeros(length(freqs), 6);
B_interp_all = zeros(nFreqInterp, 6);
K_all = zeros(nTimeSteps,6);        % Retardation functions

for DOF = 1:6
    A_eq(DOF) = vessel.powerBased.A_eq(DOF,DOF);
    B_eq(DOF) = vessel.powerBased.B_eq(DOF,DOF);
    Bv(DOF) = vessel.Bv(DOF,DOF,1);
    
    A_w = squeeze(vessel.A(DOF,DOF,:,1));
    B_w = squeeze(vessel.B(DOF,DOF,:,1));
    B_interp = interp1(freqs, B_w, freqs_uniform, 'pchip','extrap');
    B_inf = B_interp(end);

    A_w_all(:,DOF) = A_w;
    B_w_all(:,DOF) = B_w;
    B_interp_all(:,DOF) = B_interp;

    % Compute Memory Kernel K(t)
    K = zeros(nTimeSteps, 1);
    df = freqs_uniform(2) - freqs_uniform(1);
    for k = 1:nTimeSteps
        K(k) = (2/pi) * sum((B_interp-B_inf) .* cos(freqs_uniform * t(k))) * df;
        if t(k) > 50
            K(k) = 0;
        end
    end
    K_all(:,DOF) = K;  % Store K(t)

    % Cummins Equation
    M = vessel.MRB(DOF,DOF) + vessel.A(DOF,DOF,end);
    C = vessel.C(DOF,DOF);
    F_ext = waveData(:,DOF);
    eta_dof = zeros(nTimeSteps,1); % Temporary DOF result

    for k = 2:nTimeSteps-1
        tau = t(1:k);
        dtau = t(k) - tau;
        K_interp = interp1(t, K, dtau, 'linear', 0);
        memory_effect = trapz(tau(:), (K_interp(:) .* eta_dot(1:k,DOF)));
        eta_ddot(k,DOF) = (F_ext(k) -  C * eta_dof(k) - memory_effect ...
            - (B_inf + Bv(DOF)) * eta_dot(k,DOF)) / M;
        eta_dot(k+1,DOF) = eta_dot(k,DOF) + h * eta_ddot(k,DOF);
        eta_dof(k+1) = eta_dof(k) + h * eta_dot(k+1,DOF);
    end
    eta_cummins(:,DOF) = eta_dof(1:length(t));

    % Maneuvering approximation using A_eq and B_eq
    M_eq = vessel.MRB(DOF,DOF) + A_eq(DOF);
    Bv(DOF) = vessel.Bv(DOF,DOF,1);
    A_sys = [0 1;
        -C/M_eq  -(B_eq(DOF)+Bv(DOF))/M_eq];
    B_sys = [0;
        1/M_eq];
    C_sys = [1 0];
    D_sys = 0;

    % Exact ZOH discretization
    Aug = [A_sys B_sys;
        0     0     0];
    Phi = expm(Aug*h);

    Ad = Phi(1:2,1:2);
    Bd = Phi(1:2,3);

    % Initial [position; velocity] for this DOF
    xk = [0, 0]';
    eta_eq(1,DOF) = xk(1);
    eta_dot_eq(1,DOF) = xk(2);

    for k = 2:nTimeSteps-1
        xk = Ad*xk + Bd*F_ext(k);
        eta_eq(k+1,DOF) = xk(1);
        eta_dot_eq(k+1,DOF) = xk(2);
    end

end

% Convert angles to degrees for DOFs 4–6
eta_cummins(:,4:6) = rad2deg(eta_cummins(:,4:6));
eta_eq(:,4:6) = rad2deg(eta_eq(:,4:6));
eta_dot(:,4:6) = rad2deg(eta_dot(:,4:6));
eta_dot_eq(:,4:6) = rad2deg(eta_dot_eq(:,4:6));

%% Plot Results
figure(1);
dofNames = {'wave force in surge (N)', 'wave force in sway (N)', 
    'wave force in heave (N)', 'wave moment in roll (Nm)', 
    'wave moment in pitch (Nm)', 'wave moment in yaw (Nm)'};

for DOF = 1:6
    % Left column: Retardation function
    subplot(6,2,2*DOF - 1)
    plot(t, K_all(:,DOF), 'b', 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel(['K_{' num2str(DOF) num2str(DOF) '}']);
    title(['Retardation function  K_{' num2str(DOF) num2str(DOF) '}']);
    grid on;

    % Right column: 1st-order wave force
    subplot(6,2,2*DOF)
    plot(t, waveData(:,DOF), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel(['\tau_{' num2str(DOF) '}']);
    title(['1st-order ' dofNames{DOF}]);
    grid on;
end

set(findall(gcf,'type','text'),'FontSize',11)
set(findall(gcf,'type','legend'),'FontSize',8)

figure(2);
for DOF = 1:6
    % --- Added Mass subplot (left column) ---
    subplot(6,2,2*DOF - 1)
    plot(freqs, A_w_all(:,DOF), 'rx', ...
         freqs_uniform, A_eq(DOF)*ones(length(freqs_uniform),1), 'b', 'LineWidth', 1.5)
    title(['Added mass A_{' num2str(DOF) num2str(DOF) '}(ω)']);
    xlabel('Frequency (rad/s)')
    legend('A(ω)', 'A_{eq}', 'Location', 'best');
    grid on;

    % --- Damping subplot (right column) ---
    subplot(6,2,2*DOF)
    plot(freqs, B_w_all(:,DOF), 'rx', ...
         freqs_uniform, B_eq(DOF)*ones(length(freqs_uniform),1), 'b', 'LineWidth', 1.5)
    title(['Potential damping B_{' num2str(DOF) num2str(DOF) '}(ω)']);
    xlabel('Frequency (rad/s)')
    legend('B(ω)', 'B_{eq}', 'Location', 'best');
    grid on;
end

set(findall(gcf,'type','text'),'FontSize',11)
set(findall(gcf,'type','legend'),'FontSize',8)

figure(3);

% Add cruise speed to surge‐velocity perturbations
eta_dot(:,1)     = U + eta_dot(:,1);
eta_dot_eq(:,1)  = U + eta_dot_eq(:,1);

% Mapping of DOFs to subplot positions
velDOFs = [1 2 6];          % velocities → left column
posDOFs = [3 4 5];          % positions  → right column

%% -----  LEFT COLUMN : velocities (DOFs 1,2,6)  -----
for k = 1:3
    DOF = velDOFs(k);
    subplot(3,2,(k-1)*2 + 1);                       % (row k, col 1)

    plot(t, eta_dot(:,DOF),    'b-.', ...
         t, eta_dot_eq(:,DOF), 'r-', 'LineWidth',1.5);

    ylabel('Velocity');
    legend('Cummins equation', ...
           'A_{eq} and B_{eq} approx.');

    switch DOF
        case 1, title('Surge velocity (m/s)');
        case 2, title('Sway velocity (m/s)');
        case 6, title('Yaw velocity (deg/s)');
    end

    xlabel('Time (s)');
    grid on;
end

%% -----  RIGHT COLUMN : positions (DOFs 3,4,5)  -----
for k = 1:3
    DOF = posDOFs(k);
    subplot(3,2,k*2);                              % (row k, col 2)

    plot(t, eta_cummins(:,DOF), 'b-.', ...
         t, eta_eq(:,DOF),      'r-', 'LineWidth',1.5);

    ylabel('Amplitude');
    legend('Cummins equation', ...
           'A_{eq} and B_{eq} approx.');

    switch DOF
        case 3, title('Vertical (heave) position (m)');
        case 4, title('Roll angle (deg)');
        case 5, title('Pitch angle (deg)');
    end

    xlabel('Time (s)');
    grid on;
end

% Uniform font sizing
set(findall(gcf,'type','text'),   'FontSize',11)
set(findall(gcf,'type','legend'), 'FontSize',8)
