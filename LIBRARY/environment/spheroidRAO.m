function vessel = spheroidRAO(vessel,a,b,BG,zeta_roll,zeta_pitch,zn,maxDepth,verbose)
% spheroidRAO computes 6-DOF RAOs for a submerged prolate-spheroid-shaped AUV in
% waves for multiple incoming relative wave directions beta_r = [0:10:360] deg.
% The BODY-fixed coordinate origin is assumed to be in the CG.
%
% INPUTS:
%   vessel     - Structure (optional). Existing vessel data structure to be
%                updated with computed RAOs. Use [] to create a new one.
%   a          - Semi-major axis of the prolate spheroid [m].
%   b          - Semi-minor axis of the prolate spheroid [m].
%   BG         - Vertical distance between center of buoyancy (B) and center
%                of gravity (G) [m].
%   zeta_rollh - Relative damping ratio in roll
%   zeta_pitch - Relative damping ratio in pitch
%   zn         - Submergence depth of spheroid center below free surface [m].
%   maxDepth   - Total water depth [m].
%   verbose    - Logical flag (optional, default = false).
%                true  : plots RAOs, added mass, and damping curves.
%                false : silent operation.
%
% OUTPUTS:
%   vessel     - Updated structure containing computed Response Amplitude
%                Operators (RAOs) and hydrodynamic coefficients:
%                  vessel.forceRAO.w        - Wave frequency vector [rad/s]
%                  vessel.forceRAO.headings - Wave heading angles [deg]
%                  vessel.forceRAO.amp      - |H_i(ω,β)| amplitude matrices
%                  vessel.forceRAO.phase    - Phase matrices [rad]
%                  vessel.forceRAO.Re       - Real(H_i(ω,β))
%                  vessel.forceRAO.Im       - Imag(H_i(ω,β))
%
% EXAMPLES:
%   vessel = spheroidRAO(vessel,a,b,BG,zeta_roll,zeta_pitch,zn,maxDepth);
%   vessel = spheroidRAO([],2,1,0.05,0.1,0.3,5,50,true);
%
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:       T.I. Fossen
% Data:         2025-11-12
% Revisions:

clc; close all;

if nargin < 9
    verbose = false;  % Default: no plotting
end

g   = 9.81; % Acceleration of gravity (m/s2)
rho = 1025; % Density of water (kg/m3)
depth = zn; % Submergence of spheroid [m]
h_w = maxDepth; % Total water depth [m]  

if depth > h_w
    error('The vehicle submergence must satisfy: depth < h_w where h_w is the max depth')
end

% ------------------------------------------------------------------------------
% Hydrodynamic added mass (Lamb, 1932)
% ------------------------------------------------------------------------------
r44 = 0.4;
A = imlay61(a,b,zeros(6,1),r44);
MRB = spheroid(a,b,[0; 0; 0;],[0; 0; 0]);
M = MRB + A;
m = MRB(1,1);
W = m * g;

% ------------------------------------------------------------------------------
% Hydrostatic restoring (roll and pitch)
% ------------------------------------------------------------------------------
C44 = W * BG;  
C55 = W * BG;
C = diag([0 0 0 C44 C55 0]);

% ------------------------------------------------------------------------------
% Wave setup  (finite-depth compatible)
% ------------------------------------------------------------------------------
omega = linspace(0.05,4,100);
k = dispersion_relation(omega, h_w, g);  % solve ω² = g k tanh(k h)

% Exponential attenuation of dynamic pressure with submergence
% Deep water:  exp(-k*depth)
% Finite depth: cosh(k*(h_w - depth)) / cosh(k*h_w)
decay = cosh(k.*(h_w - depth)) ./ cosh(k.*h_w);

% Base excitation amplitudes (Froude–Krylov forces)
F1 = rho * g * pi * b^2 .* decay;     % surge
F3 = rho * g * pi * b^2 .* decay;   % heave
F5 = 0.1 * F3;                    % pitch (small)

beta_list = 0:10:180;
beta_rad  = deg2rad(beta_list);

% ------------------------------------------------------------------------------
% Loop over headings and compute RAOs (A(ω), B(ω) frequency dependent)
% ------------------------------------------------------------------------------
n = length(omega);
labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
H_dir = zeros(6,n,length(beta_rad));

zeta = [NaN NaN NaN zeta_roll zeta_pitch NaN];  % surge sway heave roll pitch yaw

A_freq = zeros(6,n);
B_freq = zeros(6,n);

for ib = 1:length(beta_rad)
    br = beta_rad(ib);

    % Rotate excitation for wave heading βr
    F1p = F1 .* cos(br);
    F2p = F1 .* sin(br);
    F3p = F3;
    F5p = F5 .* cos(br);
    F4p = F5 .* sin(br);
    F = [F1p; F2p; F3p; F4p; F5p; zeros(size(F1))];

    H = zeros(6,n);

    for i = 1:6
        Mii  = M(i,i);
        Aii0 = abs(A(i,i));
        Cii  = C(i,i);

        if Cii == 0
            % ------------------------------------------------------------------
            % Non-oscillatory DOFs: surge, sway, yaw → constant A and B
            % ------------------------------------------------------------------
            Aii_omega = repmat(Aii0, 1, n);
            Bii_omega = repmat(0.05 * Mii, 1, n);  % small linear damping

        else
            % ------------------------------------------------------------------
            % Oscillatory DOFs: heave, roll, pitch → freq-dependent A, B
            % ------------------------------------------------------------------
            omega_n = sqrt(Cii / (Mii + Aii0));

            % Added-mass variation (slight decrease with frequency)
            Aii_omega = Aii0 .* (1 - 0.1 * exp(-0.5 * (omega / omega_n).^2));

            % Radiation damping bell curve around resonance
            B0 = 2 * zeta(i) * (Mii + Aii0) * omega_n;
            Bii_omega = B0 .* (omega / omega_n) .* ...
                        exp(-0.5 * ((omega - omega_n) / (1.2 * omega_n)).^2);
        end

        % Store for plotting
        A_freq(i,:) = Aii_omega;
        B_freq(i,:) = Bii_omega;

        % Compute complex RAO
        H(i,:) = F(i,:) ./ (Cii - omega.^2 .* (Mii + Aii_omega) + ...
            1i .* omega .* Bii_omega);
    end

    H_dir(:,:,ib) = H;
end

if verbose
    fprintf('\nNatural Frequencies and Damping Ratios (Approx.):\n');
    fprintf('--------------------------------------------------\n');
    for i = 1:6
        Mii = M(i,i);
        Aii0 = abs(A(i,i));
        Cii = C(i,i);
        if Cii > 0
            omega_n = sqrt(Cii / (Mii + Aii0));
            Bii1 = 2*zeta(i)*(Mii + Aii0); % from definition
            zeta_eff = Bii1 / (2*(Mii + Aii0)*omega_n);
            fprintf('%-6s: ω_n = %5.3f rad/s, ζ_eff = %.3f\n', ...
                labels{i}, omega_n, zeta_eff);
        end
    end
end

% ------------------------------------------------------------------------------
% Post-process RAOs: amplitude, phase, and complex form (0–360°)
% ------------------------------------------------------------------------------
vessel.forceRAO = struct();
vessel.forceRAO.w = omega;  % [rad/s]

numDOF = 6;

amp   = cell(1,numDOF);
phase = cell(1,numDOF);
ReH   = cell(1,numDOF);
ImH   = cell(1,numDOF);

% Define consistent 0–350° headings (36 directions)
vessel.headings = 0:10:350;  % [deg]

for i = 1:numDOF
    % Extract complex RAOs for all directions (0–180)
    H_half = squeeze(H_dir(i,:,:));   % size: n_omega × n_beta

    % Mirror to 360° using complex conjugate symmetry:
    % H(ω, β+180°) = conj(H(ω, β))
    H_full = [H_half, conj(fliplr(H_half(:,2:end-1)))];

    % Store components
    amp{i}   = abs(H_full);
    phase{i} = angle(H_full);
    ReH{i}   = real(H_full);
    ImH{i}   = imag(H_full);
end

% Store final results
vessel.forceRAO.amp   = amp;
vessel.forceRAO.phase = phase;
vessel.forceRAO.Re    = ReH;
vessel.forceRAO.Im    = ImH;

% ------------------------------------------------------------------------------
% Plot all 6 DOFs 
% ------------------------------------------------------------------------------
if verbose
    figure(1);
    for i = 1:6
        subplot(3,2,i); hold on; grid on;
        for ib = 1:19
            plot(omega, abs(H_dir(i,:,ib)), 'LineWidth', 1.1, ...
                'DisplayName', sprintf('\\beta_r = %d°', beta_list(ib)));
        end
        xlabel('\omega [rad/s]');
        ylabel(sprintf('|H_%d|', i));
        title(labels{i});
    end
    sgtitle('RAOs of a Submerged Prolate Spheroid (0-360 deg');
end

% ------------------------------------------------------------------------------
% Plot frequency-dependent Added Mass A(ω)
% ------------------------------------------------------------------------------
if verbose
    figure(2);
    titles = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};

    for i = 1:6
        subplot(3,2,i); hold on; grid on;
        plot(omega, A_freq(i,:), 'b', 'LineWidth',1.3);
        xlabel('\omega [rad/s]');
        ylabel('A_{ii}  [kg or kg·m²]');
        title(titles{i});
    end
    sgtitle('Frequency-Dependent Added Mass A(\omega)');
end

% ------------------------------------------------------------------------------
% Plot frequency-dependent Damping B(ω)
% ------------------------------------------------------------------------------
if verbose
    figure(3);
    for i = 1:6
        subplot(3,2,i); hold on; grid on;
        plot(omega, B_freq(i,:), 'b', 'LineWidth',1.3);
        xlabel('\omega [rad/s]');
        ylabel('B_{ii}  [N·s/m or N·m·s/rad]');
        title(titles{i});
    end
    sgtitle('Frequency-Dependent Damping B(\omega)');
end


% ------------------------------------------------------------------------------
% Helper function
% ------------------------------------------------------------------------------
    function k = dispersion_relation(omega, h, g)
        % Dispersion_relation  Solve ω² = g k tanh(k h) for finite water depth.
        % Uses Newton–Raphson iteration with deep-water initial guess.

        k = omega.^2 ./ g; % Initial guess
        for iter = 1:20
            f  = g.*k.*tanh(k.*h) - omega.^2;
            df = g.*tanh(k.*h) + g.*h.*k.*sech(k.*h).^2;
            k  = k - f ./ df;
        end
    end

end