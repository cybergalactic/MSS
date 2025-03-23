% This script simulates the response of a vessel in waves using Cummins' 
% equation and an equivalent maneuvering model. It calculates the wave-induced 
% forces, solves the full hydrodynamic model including memory effects, and 
% compares it with a simplified approximation using added mass and damping 
% coefficients (Aeq, Beq) according to: 
%
%   A_eq(i,j) = ∫ A(i,j,ω) S(ω) dω  /  ∫ S(ω) dω
%   B_eq(i,j) = ∫ B(i,j,ω) S(ω) dω  /  ∫ S(ω) dω
%
% where:
%   - A(i,j,ω) and B(i,j,ω) are the added mass and damping coefficients 
%     at frequency ω for each matrix element (i,j).
%   - S(ω) is the wave energy spectrum, computed using waveSpectrum().
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
rng(1); % Set random generator seed to 1 when generating stochastic waves

%% USER INPUTS
h  = 0.05; % Sampling time [s]
T_final = 200; % Final simulation time [s]
plotFlag = 0; % Set to 1 to plot 6x6 matrix elememts, 0 for no plot

load s175; % Load vessel structure
U = 0; % Speed (m/s)
psi = 0; % Heading angle (rad)
beta_wave = deg2rad(50); % Wave direction relative bow, 0 for following sea, 180 for head sea
maxFreq = 3.0; % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 60; % Number of wave frequency intervals (>50)
   
% Calculate the wave spectrum power intensity S(Omega) for each frequency
spectrumNo = 7; % JONSWAP
Hs = 3; % Significant wave height (m)
w0 = 1.2;  % Peak frequency (rad/s)
gamma = 3.3; % Peakedness factor
Parameter = [Hs, w0, gamma]; % Spectrum parameters

% Time vector from 0 to T_final     
t = 0:h:T_final;      
nTimeSteps = length(t);

% Wave spectrum, one direction
omegaMax = vessel.forceRAO.w(end);  % Max frequency in RAO dataset

[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumNo, ...
    Parameter, numFreqIntervals, omegaMax);

% 6-DOF generalized wave forces using firts-order force RAOs
waveData = zeros(nTimeSteps,7); % Pre-allocate table
for i = 1:nTimeSteps
    [tau_wave1, waveElevation] = waveForceRAO(t(i), ...
        S_M, Amp, Omega, mu, vessel, U, psi, beta_wave, numFreqIntervals);
    waveData(i,:) = [tau_wave1' waveElevation];
end

%% Compute Aeq and Beq for heave, roll and pitch
vessel.B = vessel.B + vessel.Bv; % Add viscous damping

for DOF = 3:5
    vessel = computeManeuveringModel(vessel, 1, w0, plotFlag);
    A_eq = vessel.A_eq(DOF,DOF);
    B_eq = vessel.B_eq(DOF,DOF);

    %% Compute Memory Function K(t) for Cummins Equation
    A_w = squeeze(vessel.A(DOF,DOF,:,1));

    % Interpolate B(w) to Evenly Spaced Frequency Grid
    freqs = vessel.freqs;
    nFreqInterp = 200; % Higher resolution
    freqs_uniform = linspace(min(freqs), max(freqs), nFreqInterp)';
    B_w = squeeze(vessel.B(DOF,DOF,:,1));
    B_interp = interp1(freqs, B_w, freqs_uniform, 'pchip','extrap');
    B_inf = B_interp(end);

    %% Compute Memory Function K(t)
    K = zeros(nTimeSteps, 1); % Initialize K(t)
    df = freqs_uniform(2) - freqs_uniform(1);

    for k = 1:nTimeSteps
        K(k) = (2/pi) * sum((B_interp-B_inf) .* cos(freqs_uniform * t(k))) * df;
        if t(k) > 50
            K(k) = 0;
        end
    end

    %% Solve Cummins Equation (Full Model)
    M = vessel.MRB(DOF,DOF) + vessel.A(DOF,DOF,end); % Total inertia (includes A_inf)
    C = vessel.C(DOF,DOF);

    % External excitation force, 1st-order wave loads 
    F_ext = waveData(:,DOF);

    eta_cummins = zeros(nTimeSteps,1);  % Displacement
    eta_dot = zeros(nTimeSteps,1);      % Velocity
    eta_ddot = zeros(nTimeSteps,1);     % Acceleration

    for k = 2:nTimeSteps
        % Compute the time differences
        tau = t(1:k);         % Past time values
        dtau = t(k) - tau;    % Time differences (t - tau)

        % Interpolate K at the required (t - tau) values
        K_interp = interp1(t, K, dtau, 'linear', 0); 

        % Compute the convolution integral for memory effect
        memory_effect = trapz(tau(:), (K_interp(:) .* eta_dot(1:k)));

        % Compute acceleration (Newton’s Second Law) at time k
        eta_ddot(k) = (F_ext(k) - C * eta_cummins(k) - memory_effect ...
            - B_inf * eta_dot(k)) / M;

        % Integrate velocity and displacement using Euler's method
        eta_dot(k+1) = eta_dot(k) + h * eta_ddot(k);   % Velocity update
        eta_cummins(k+1) = eta_cummins(k) + h * eta_dot(k+1); % Displacement update
    end

    eta_cummins = eta_cummins(1:length(t));  % Trim to match the length of t

    %% Solve Using Aeq and Beq (Maneuvering Model)
    M_eq = vessel.MRB(DOF,DOF) + A_eq;  % Effective mass for simplified model

    % Define state-space system: [dx/dt] = Ax + Bu, y = Cx + Du
    A_sys = [0 1; -C/M_eq -B_eq/M_eq]; % State matrix
    B_sys = [0; 1/M_eq];               % Input matrix
    C_sys = [1 0];                     % Output matrix (displacement)
    D_sys = 0;                         % Direct transmission

    sys_eq = ss(A_sys, B_sys, C_sys, D_sys); % Create system

    % Simulate response
    eta_eq = lsim(sys_eq, F_ext, t);

    if DOF == 4 || DOF == 5 
        eta_cummins = rad2deg(eta_cummins);
        eta_eq = rad2deg(eta_eq);
    end

    %% Plot Results
    DOFtext1 = {'Heave Position (m)', 'Roll Angle (deg)', 'Pitch Angle (deg)'};
    DOFtext2 = {'Heave Retardation Function K_{33}(t)', ...
        'Roll Retardation Function K_{44}(t)', 'Pitch Retardation Function K_{55}(t)'};
    DOFtext3 = {'A_{33}(ω)', 'A_{44}(ω)', 'A_{55}(ω)'};
    DOFtext4 = {'B_{33}(ω)', 'B_{44}(ω)', 'B_{55}(ω)'};

    figure(1);
    subplot(3,1,DOF-2)
    plot(t, eta_cummins, 'k', t, eta_eq, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Cummins Equation', 'Aeq and Beq Approximation');
    title(DOFtext1{DOF-2});
    grid;
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    set(findall(gcf,'type','line'),'linewidth',2)

    figure(2);
    subplot(3,1,DOF-2)
    plot(t,K, 'LineWidth', 2)
    xlabel('Time (s)');
    ylabel('Memory Kernel K(t)');
    title(DOFtext2{DOF-2});
    grid;
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','line'),'linewidth',2)

    figure(3);
    subplot(3,1,DOF-2)
    plot(freqs, A_w,'rx', ...
        freqs_uniform,A_eq*ones(length(freqs_uniform),1),'b','LineWidth', 2)
    title(DOFtext3{DOF-2});
    legend('A(ω)','A_{eq}');
    grid;
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    set(findall(gcf,'type','line'),'linewidth',2)

    figure(4);
    subplot(3,1,DOF-2)
    plot(freqs_uniform, B_interp,'g', ...
        freqs, B_w,'rx',...
        freqs_uniform, B_eq*ones(length(freqs_uniform),1), 'b', ...
        'LineWidth', 2)
    title(DOFtext4{DOF-2});
    legend('Interpolated','B(ω)','B_{eq}');
    grid;
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    set(findall(gcf,'type','line'),'linewidth',2)

end