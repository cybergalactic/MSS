% This script simulates the response of a vessel in waves using Cummins' 
% equation and an equivalent maneuvering model. It calculates the wave-induced 
% forces, solves the full hydrodynamic model including memory effects, and 
% compares it with a simplified approximation using equivalent added mass and 
% damping coefficients (A_eq, B_eq) according to: 
%
%   A_eq(i,j,ω_p,velocity) = ∫ A(i,j,ω,velocity) S(ω,ω_p(k)) dω / ∫ S(ω,ω_p(k)) dω
%   B_eq(i,j,ω_p,velocity) = ∫ B(i,j,ω,velocity) S(ω,ω_p(k)) dω / ∫ S(ω,ω_p(k)) dω
%
% where:
%   - A(i,j,ω,velocity) and B(i,j,ω,velocity) are the added mass and damping 
%     coefficients at frequency and velocity for each matrix element (i,j).
%   - S(ω,ω_p) is the wave energy spectrum.
%   - The integrals are evaluated numerically using trapezoidal integration.
%
% This ensures that the kinetic energy and power dissipation properties 
% of the frequency-dependent system are preserved in the equivalent 
% constant matrices.
%
% Author:    Thor I. Fossen
% Date:      2025-03-10
% Revisions: 

clear waveForceRAO; % Clear persistent RAO tables
clearvars; 
close all;
rng(1); % Set random generator seed to 1 when generating stochastic waves

%% SIMULATOR CONFIGURATION
h  = 0.05; % Sampling time (s)
T_final = 200; % Final simulation time (s)
plotFlag = 0; % Set to 1 to plot 6x6 matrix elememts, 0 for no plot
vesselChoice = 3; % Choose vessel type 1, 2, 3

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
beta_wave = deg2rad(50); % Wave direction, 0 for following sea, 180 for head sea
maxFreq = 3.0; % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 60; % Number of wave frequency intervals (>50)
   
% Calculate the wave spectrum power intensity S(Omega) for each frequency
spectrumType = 'JONSWAP'; 
Hs = 4; % Significant wave height (m)
omega_p = 0.8;  % Wave spectrum peak frequencies (rad/s)
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

%% Compute Aeq and Beq 
vessel = computeManeuveringModel(vessel, omega_p, plotFlag);

%% Compute Cummins and Maneuvering Model Responses
freqs = vessel.freqs;
nFreqInterp = 200;
freqs_uniform = linspace(min(freqs), max(freqs), nFreqInterp)';
vessel.B = vessel.B + vessel.Bv; % Add viscous damping

% Initialize storage for all DOFs
eta_cummins = zeros(nTimeSteps,6);  % Displacement (Cummins)
eta_eq = zeros(nTimeSteps,6);       % Displacement (Aeq-Beq)
eta_dot = zeros(nTimeSteps,6);      % Velocity (Cummins)
eta_ddot = zeros(nTimeSteps,6);     % Acceleration (Cummins)
eta_dot_eq = zeros(nTimeSteps,6);   % Velocity (Maneuvering)
A_eq = zeros(6,1);  
B_eq = zeros(6,1);

A_w_all = zeros(length(freqs), 6);
B_w_all = zeros(length(freqs), 6);
B_interp_all = zeros(nFreqInterp, 6);
K_all = zeros(nTimeSteps,6);        % Retardation functions

for DOF = 1:6
    A_eq(DOF) = vessel.A_eq(DOF,DOF);
    B_eq(DOF) = vessel.B_eq(DOF,DOF);

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
        eta_ddot(k,DOF) = (F_ext(k) - C * eta_dof(k) - memory_effect ...
            - B_inf * eta_dot(k,DOF)) / M;
        eta_dot(k+1,DOF) = eta_dot(k,DOF) + h * eta_ddot(k,DOF);
        eta_dof(k+1) = eta_dof(k) + h * eta_dot(k+1,DOF);
    end
    eta_cummins(:,DOF) = eta_dof(1:length(t));

    % Maneuvering Approximation using A_eq and B_eq
    M_eq = vessel.MRB(DOF,DOF) + A_eq(DOF);
    A_sys = [0 1; -C/M_eq -B_eq(DOF)/M_eq];
    B_sys = [0; 1/M_eq];
    C_sys = [1 0];
    D_sys = 0;
    sys_eq = ss(A_sys, B_sys, C_sys, D_sys);
    [eta_eq(:,DOF), ~, x_eq_states]  = lsim(sys_eq, F_ext, t); % Position
    eta_dot_eq(:,DOF) = x_eq_states(:,2); % Velocity 
end

% Convert angles to degrees for DOFs 4–6
eta_cummins(:,4:6) = rad2deg(eta_cummins(:,4:6));
eta_eq(:,4:6) = rad2deg(eta_eq(:,4:6));
eta_dot(:,4:6) = rad2deg(eta_dot(:,4:6));
eta_dot_eq(:,4:6) = rad2deg(eta_dot_eq(:,4:6));

%% Plot Results
figure(1);
dofNames = {'Wave Force in Surge (N)', 'Wave Force in Sway (N)', 
    'Wave Force in Heave (N)', 'Wave Moment in Roll (Nm)', 
    'Wave Moment in Pitch (Nm)', 'Wave Moment in Yaw (Nm)'};

for DOF = 1:6
    % Left column: Retardation function
    subplot(6,2,2*DOF - 1)
    plot(t, K_all(:,DOF), 'b', 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel(['K_{' num2str(DOF) num2str(DOF) '}']);
    title(['Retardation Function  K_{' num2str(DOF) num2str(DOF) '}']);
    grid on;

    % Right column: 1st-order wave force
    subplot(6,2,2*DOF)
    plot(t, waveData(:,DOF), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel(['\tau_{' num2str(DOF) '}']);
    title(['1st-Order ' dofNames{DOF}]);
    grid on;
end

set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',10)

figure(2);
for DOF = 1:6
    % --- Added Mass subplot (left column) ---
    subplot(6,2,2*DOF - 1)
    plot(freqs, A_w_all(:,DOF), 'rx', ...
         freqs_uniform, A_eq(DOF)*ones(length(freqs_uniform),1), 'b', 'LineWidth', 1.5)
    title(['Added Mass A_{' num2str(DOF) num2str(DOF) '}(ω)']);
    legend('A(ω)', 'A_{eq}', 'Location', 'best');
    grid on;

    % --- Damping subplot (right column) ---
    subplot(6,2,2*DOF)
    plot(freqs_uniform, B_interp_all(:,DOF), 'g', ...
         freqs, B_w_all(:,DOF), 'rx', ...
         freqs_uniform, B_eq(DOF)*ones(length(freqs_uniform),1), 'b', 'LineWidth', 1.5)
    title(['Damping B_{' num2str(DOF) num2str(DOF) '}(ω)']);
    legend('Interpolated', 'B(ω)', 'B_{eq}', 'Location', 'best');
    grid on;
end

set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',10)

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

    plot(t, eta_dot(:,DOF),    'b', ...
         t, eta_dot_eq(:,DOF), 'r', 'LineWidth',1.5);

    ylabel('Velocity');
    legend('Cummins Equation', ...
           'A_{eq} and B_{eq} approx.');

    switch DOF
        case 1, title('Surge Velocity (m/s)');
        case 2, title('Sway Velocity (m/s)');
        case 6, title('Yaw Velocity (deg/s)');
    end

    xlabel('Time (s)');
    grid on;
end

%% -----  RIGHT COLUMN : positions (DOFs 3,4,5)  -----
for k = 1:3
    DOF = posDOFs(k);
    subplot(3,2,k*2);                              % (row k, col 2)

    plot(t, eta_cummins(:,DOF), 'b', ...
         t, eta_eq(:,DOF),      'r', 'LineWidth',1.5);

    ylabel('Amplitude');
    legend('Cummins Equation', ...
           'A_{eq} and B_{eq} approx.');

    switch DOF
        case 3, title('Vertical (Heave) Position (m)');
        case 4, title('Roll Angle (deg)');
        case 5, title('Pitch Angle (deg)');
    end

    xlabel('Time (s)');
    grid on;
end

% Uniform font sizing
set(findall(gcf,'type','text'),   'FontSize',12)
set(findall(gcf,'type','legend'), 'FontSize',10)
