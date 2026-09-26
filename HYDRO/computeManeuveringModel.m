function vessel = computeManeuveringModel(vessel,  omega_p, ...
    kappa_126, delta_zeta_345, plotFlag)
% computeManeuveringModel is compatible with MATLAB and GNU Octave (www.octave.org).
% Computes the power-based equivalent added mass A_eq and potential damping B_eq 
% by integrating the frequency-dependent hydrodynamic matrices A_U(omega) and 
% B_U(omega) using the wave spectrum S(omega)as a weighting function. The 
% diagonal viscous damping matrix Bv is computed from damping increments 
% for surge, sway, and yaw and damping-ratio increments for heave, roll, and 
% pitch. The total effective damping matrix is D = B_eq + Bv.
%
% The equivalent matrices are computed as:
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
% Inputs:
%   vessel         - Structure containing vessel hydrodynamic data
%   omega_p        - Wave peak frequency (rad/s)
%   kappa_126      - Relative viscous damping increments for DOFs 1, 2 and 6
%                    (default: [0.05 0.05 0.05])
%   delta_zeta_345 - Viscous damping-ratio increments for DOFs 3, 4 and 5
%                    (default: [0 0.1 0])
%   plotFlag       - Set to 1 to plot A(omega) and B(omega), 0 otherwise
%
% Outputs:
%   vessel.powerBased.omega_p - Wave spectrum peak frequency
%   vessel.powerBased.eq.A_eq - Equivalent added mass matrix
%   vessel.powerBased.B_eq    - Equivalent damping matrix
%
% Example call:
%   load supply; % Other vessels: s175, tanker, fpso, semisub,
%                                 capytaineTestShip
%
%   omega_p = 1.0
%
%   vessel = computeManeuveringModel(vessel, omega_p);
%   vessel = computeManeuveringModel(vessel, omega_p, ...
%       kappa_126, delta_zeta_345);
%   vessel = computeManeuveringModel(vessel, omega_p, ...
%       kappa_126, delta_zeta_345, 1);
%   vessel = computeManeuveringModel(vessel, omega_p, [], [], 1);
%
%   disp(vessel.powerBased.A_eq)
%   disp(vessel.powerBased.B_eq)
%
% Author: Thor I. Fossen
% Date: 2025-03-10
% Revisions: 
%   2026-04-07 Use only potential damping vessel.B when computing B_eq.
%   2026-09-26 Introduced the structure vessel.powerBased and added 
%      formulas for the viscous damping matrix Bv.

% Default damping increments (DOFs 1-2-6)
if nargin < 3 || isempty(kappa_126)
    kappa_126 = [0.05 0.05 0.05];  % 5 percent increase of damping
end

% Default damping-ratio increments (DOFs 3-4-5)
if nargin < 4 || isempty(delta_zeta_345)
    delta_zeta_345 = [0 0.1 0]; % Increase the damping factor in roll by 0.1
end

% Default plot flag
if nargin < 5 || isempty(plotFlag)
    plotFlag = 0;
end

% Check number of velocity cases
if isfield(vessel, 'velocities') && ~isempty(vessel.velocities)
    nvel = length(vessel.velocities);
else
    nvel = 1;
end

A_all = vessel.A; % Added mass
B_all = vessel.B; % Potential damping 

% Frequency data used for power-based averaging
freqs = vessel.freqs;

% Exclude artificial frequency omega = 10 rad/s representing infinity
idx = freqs < 10;
freqs = freqs(idx);

omega_min = min(freqs);
omega_max = max(freqs);

% Avoid omega = 0 to prevent numerical issues in spectrum normalization
if omega_min == 0
    omega_min = 1e-6;
end

% Define finer frequency grid for interpolation
freqs_fine = linspace(omega_min, omega_max, 100)';
nOmega = length(omega_p);

% PM wave spectrum parameters
alpha = 8.1e-3 * (9.81)^2;
beta = 0.74;

% Initialize equivalent matrices of dimension [6, 6, nOmega, nvel]
Aeq_all = zeros(6, 6, nOmega, nvel);
Beq_all = zeros(6, 6, nOmega, nvel);

% Loop over all velocities
for velNo = 1:nvel
    A_w = A_all(:,:,idx,velNo);
    B_w = B_all(:,:,idx,velNo);

    % Zero spectral moment m_0
    S = alpha ./ freqs_fine.^5 .* exp(-beta * (omega_p ./ freqs_fine).^4);
    m_0 = trapz(freqs_fine, S);

    % Normalized wave spectrum
    S_N = S / m_0;

    % Loop over DOFs
    for i = 1:6
        for j = 1:6
            A_ij_w = squeeze(A_w(i,j,:));
            B_ij_w = squeeze(B_w(i,j,:));

            A_interp = interp1(freqs, A_ij_w, freqs_fine, 'pchip');
            B_interp = interp1(freqs, B_ij_w, freqs_fine, 'pchip');

            Aeq_all(i,j,velNo) = trapz(freqs_fine, A_interp .* S_N);
            Beq_all(i,j,velNo) = trapz(freqs_fine, B_interp .* S_N);
        end
    end

end

% Store in vessel
vessel.powerBased.omega_p = omega_p;
vessel.powerBased.A_eq = Aeq_all;
vessel.powerBased.B_eq = Beq_all;


%% Power-based model matrices
vessel.powerBased.Bv = zeros(6);
vessel.powerBased.G  = zeros(6);

% Total inertia matrix
vessel.powerBased.M = vessel.MRB + vessel.powerBased.A_eq;

% Restoring matrix for heave, roll and pitch
vessel.powerBased.G([3 4 5],[3 4 5]) = ...
    vessel.C([3 4 5],[3 4 5],1);

% DOFs 1, 2 and 6: relative viscous damping increments
idx = [1 2 6];
for k = 1:3
    i = idx(k);
    vessel.powerBased.Bv(i,i) = ...
        kappa_126(k) * vessel.powerBased.B_eq(i,i);
end

% DOFs 3, 4 and 5: viscous damping-ratio increments
idx = [3 4 5];
for k = 1:3
    i = idx(k);
    vessel.powerBased.Bv(i,i) = ...
        2 * delta_zeta_345(k) * ...
        sqrt(vessel.powerBased.M(i,i) * vessel.powerBased.G(i,i));
end

%% Store final results
vessel.powerBased.D = vessel.powerBased.B_eq + vessel.powerBased.Bv;
vessel.powerBased.MA  = vessel.powerBased.A_eq;
vessel.powerBased.MRB = vessel.MRB;

vessel.powerBased.kappa_126 = kappa_126;
vessel.powerBased.delta_zeta_345 = delta_zeta_345;

vessel.powerBased.T_126 = zeros(1,3);
idx = [1 2 6];

for k = 1:3
    i = idx(k);
    vessel.powerBased.T_126(k) = ...
        vessel.powerBased.M(i,i) / vessel.powerBased.D(i,i);
end

%% Optional plotting for velocity #1
if plotFlag == 1
    plotAB_eq(vessel, 'A', 1);
    plotAB_eq(vessel, 'B', 1);
end

end