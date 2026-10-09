function vessel = computeManeuveringModel(vessel, omega_p, ...
    aperiodicDamping, delta_zeta, plotFlag)
% Computes the power-based equivalent added mass A_eq and potential damping B_eq
% by integrating the frequency-dependent hydrodynamic matrices A_U(omega) and
% B_U(omega) using the wave spectrum S(omega) as a weighting function. The
% diagonal viscous damping matrix Bv is computed using one of two call patterns:
%
%   Floating vessel:   kappa_126        (3 elements) and
%                      delta_zeta_345   (3 elements)
%   Submerged vehicle: T_1236           (4 elements) and
%                      delta_zeta_45    (2 elements)
%
% Floating-vessel formulation:
%   kappa_126 = [kappa_1 kappa_2 kappa_6] contains dimensionless viscous-
%   damping increments for the unrestrained surge, sway and yaw modes:
%
%       Bv(i,i) = kappa_i * B_eq(i,i),       i = 1, 2, 6.
%
%   Hence, kappa_i = 0.05 adds viscous damping equal to 5 percent of the
%   equivalent potential damping in DOF i. The vector
%   delta_zeta_345 = [delta_zeta_3 delta_zeta_4 delta_zeta_5] specifies
%   increments in the heave, roll and pitch damping ratios:
%
%       Bv(i,i) = 2 * delta_zeta_i * sqrt(M(i,i) * G(i,i)), i = 3, 4, 5.
%
% Submerged-vehicle formulation:
%   T_1236 = [T_1 T_2 T_3 T_6] specifies positive target time constants
%   [s] for the unrestrained surge, sway, heave and yaw modes. For the
%   uncoupled scalar model M(i,i)*nu_dot_i + D(i,i)*nu_i = 0,
%
%       T_i = M(i,i) / D(i,i),
%       Bv(i,i) = M(i,i) / T_i - B_eq(i,i),  i = 1, 2, 3, 6.
%
%   Thus, smaller T_i gives greater damping. The requested T_i must not
%   require negative Bv(i,i). In a coupled model, the modal time constants
%   can differ slightly from these diagonal, uncoupled target values. The
%   vector delta_zeta_45 = [delta_zeta_4 delta_zeta_5] specifies damping-
%   ratio increments for the restored roll and pitch modes using the same
%   formula as above.
%
% The vector dimensions select the formulation (3+3 or 4+2 parameters).
% The selection is then checked using the heave restoring coefficient
% G(3,3). For a surface-piercing floating vessel, a vertical displacement
% changes the displaced volume through the waterplane area and produces a
% linear hydrostatic heave restoring force; consequently, G(3,3) is nonzero.
% For a freely submerged, constant-volume vehicle in homogeneous water, a
% small vertical displacement changes neither weight nor buoyancy, so there
% is no linear hydrostatic heave stiffness and G(3,3) = 0. Heave is therefore
% an aperiodic mode and requires T_3. Numerically, G(3,3) is treated as zero
% when abs(G(3,3)) <= 1e-10 * max(1,norm(G,'fro')).
% This classification assumes no tether, vertical spring, or modeled
% depth-dependent buoyancy that would give a submerged vehicle heave stiffness.
%
% The total effective damping matrix is D = B_eq + Bv. The corresponding system
% inertia matrix is M = MRB + MA where MA = A_eq.
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
%   omega_p        - Wave peak frequency (rad/s), typically 0.8 to 1.0 rad/s
%   aperiodicDamping - Floating vessel: kappa_126, the dimensionless relative
%                     viscous-damping increments for DOFs 1, 2 and 6. The
%                     default is [0.05 0.05 0.05].
%                     Submerged vehicle: T_1236, the target time constants
%                     [s] for DOFs 1, 2, 3 and 6. The default is
%                     [50 5 5 5] s.
%   delta_zeta       - Floating vessel: damping-ratio increments
%                      delta_zeta_345 for DOFs 3, 4 and 5. The default is
%                      [0 0.1 0].
%                      Submerged vehicle: damping-ratio increments
%                      delta_zeta_45 for DOFs 4 and 5. The default is [0.2 0.2].
%   plotFlag       - Set to 1 to plot A(omega) and B(omega), 0 otherwise
%
% Outputs:
%   vessel.powerBased.omega_p - Wave spectrum peak frequency used for
%                               power-based weighting
%   vessel.powerBased.A_eq    - 6x6x equivalent added mass matrix
%   vessel.powerBased.B_eq    - 6x6 equivalent damping matrix
%   vessel.powerBased.Bv      - 6x6 viscous damping matrix
%
% Usage:
%   load <vessel>;  % supply, s175, tanker, fpso, semisub, testShip, LAUV_marie
%
%   If powerBased inputs have not previously been stored in the vessel data
%   structure, you will be asked to supply omega_p and the parameters for
%   viscous damping:
%
%     vessel = computeManeuveringModel(vessel);
%     vessel = computeManeuveringModel(vessel, omega_p);
%     vessel = computeManeuveringModel(vessel, omega_p, [], [], 1); % Plotting
%
%   Supply all model inputs explicitly (no prompts are displayed):
%
%     vessel = computeManeuveringModel(vessel, omega_p, kappa_126, delta_zeta_345);
%     vessel = computeManeuveringModel(vessel, 0.8, [0.05 0.05 0.05], [0 0.1 0]);
%     vessel = computeManeuveringModel(vessel, omega_p, T_1236, delta_zeta_45);
%     vessel = computeManeuveringModel(vessel, 0.8, [50 5 5 5], [0.2 0.2]);
%
%   Supply all model inputs explicitly and enable plotting.
%
%     vessel = computeManeuveringModel(vessel, omega_p, kappa_126, delta_zeta_345, 1);
%
% Display the matrices:
%
%     display(vessel.powerBased.A_eq, 'A_eq')
%     display(vessel.powerBased.B_eq, 'B_eq')
%     display(vessel.powerBased.Bv, 'Bv')
%
% Reference:
%     Fossen, T. I. (2025). Maneuvering Coefficient Estimation from Frequency-
%     Dependent Added Mass and Damping: A Power-Based Approach.
%     Ocean Engineering, 341, 122494.
%
% Author: Thor I. Fossen
% Date: 2025-03-10
% Revisions: 
%   2026-04-07 Use only potential damping vessel.B when computing B_eq.
%   2026-09-27 Introduced the structure vessel.powerBased and added
%      formulas for the viscous damping matrix Bv. Prompt for damping increments
%      when they are not supplied.
%   2026-10-04 Added the submerged-vehicle formulation using T_1236 and
%      delta_zeta_45. The G(3,3) restoring coefficient validates the call.

% ------------------------------------------------------------------------------
% Wave-spectrum peak frequency. Reuse a saved value when the argument is
% omitted; an explicitly empty argument requests a new value.
% ------------------------------------------------------------------------------
if nargin < 2
    if isfield(vessel, 'powerBased') && ...
            isfield(vessel.powerBased, 'omega_p') && ...
            ~isempty(vessel.powerBased.omega_p)
        omega_p = vessel.powerBased.omega_p;
    else
        omega_p = input('Wave-spectrum peak frequency omega_p (rad/s): ');
    end
elseif isempty(omega_p)
    omega_p = input('Wave-spectrum peak frequency omega_p (rad/s): ');
end

% ------------------------------------------------------------------------------
% Restoring matrix and heave-mode classification. A floating vessel has
% hydrostatic heave restoring, whereas a fully submerged vehicle does not.
% ------------------------------------------------------------------------------
G = zeros(6);
G([3 4 5],[3 4 5]) = vessel.C([3 4 5],[3 4 5],1);
tolG = 1e-10 * max(1, norm(G, 'fro'));
heaveIsRestored = abs(G(3,3)) > tolG;

% ------------------------------------------------------------------------------
% Reuse saved damping inputs when omitted. An explicitly empty argument
% displays the default prompt for the formulation identified by G(3,3).
% ------------------------------------------------------------------------------
if nargin < 3
    if heaveIsRestored && isfield(vessel, 'powerBased') && ...
            isfield(vessel.powerBased, 'kappa_126') && ...
            ~isempty(vessel.powerBased.kappa_126)
        aperiodicDamping = vessel.powerBased.kappa_126;
    elseif ~heaveIsRestored && isfield(vessel, 'powerBased') && ...
            isfield(vessel.powerBased, 'T_1236') && ...
            ~isempty(vessel.powerBased.T_1236)
        aperiodicDamping = vessel.powerBased.T_1236;
    else
        aperiodicDamping = [];
    end
end

if isempty(aperiodicDamping)
    if heaveIsRestored
        defaultValue = [0.05 0.05 0.05];
        prompt = ['Relative viscous damping increments kappa_126 for ' ...
            'DOFs 1, 2, 6'];
    else
        defaultValue = [50 5 5 5];
        prompt = 'Aperiodic time constants T_1236 [s] for DOFs 1, 2, 3, 6';
    end
    aperiodicDamping = input(sprintf('%s (default: %s): ', ...
        prompt, mat2str(defaultValue)));
    if isempty(aperiodicDamping)
        aperiodicDamping = defaultValue;
    end
end

if nargin < 4
    if heaveIsRestored && isfield(vessel, 'powerBased') && ...
            isfield(vessel.powerBased, 'delta_zeta_345') && ...
            ~isempty(vessel.powerBased.delta_zeta_345)
        delta_zeta = vessel.powerBased.delta_zeta_345;
    elseif ~heaveIsRestored && isfield(vessel, 'powerBased') && ...
            isfield(vessel.powerBased, 'delta_zeta_45') && ...
            ~isempty(vessel.powerBased.delta_zeta_45)
        delta_zeta = vessel.powerBased.delta_zeta_45;
    else
        delta_zeta = [];
    end
end

if isempty(delta_zeta)
    if heaveIsRestored
        defaultValue = [0 0.1 0];
        prompt = ['Viscous damping-ratio increments delta_zeta_345 for ' ...
            'DOFs 3, 4, 5'];
    else
        defaultValue = [0.2 0.2];
        prompt = ['Viscous damping-ratio increments delta_zeta_45 for ' ...
            'DOFs 4, 5'];
    end
    delta_zeta = input(sprintf('%s (default: %s): ', ...
        prompt, mat2str(defaultValue)));
    if isempty(delta_zeta)
        delta_zeta = defaultValue;
    end
end

% Vector dimensions select the call pattern; G(3,3) verifies it.
surfaceCall = numel(aperiodicDamping) == 3 && numel(delta_zeta) == 3;
submergedCall = numel(aperiodicDamping) == 4 && numel(delta_zeta) == 2;
if ~surfaceCall && ~submergedCall
    error(['Use 3+3 damping parameters (kappa_126, delta_zeta_345) for a ' ...
        'floating vessel or 4+2 parameters (T_1236, delta_zeta_45) for a ' ...
        'submerged vehicle.']);
end
if surfaceCall && ~heaveIsRestored
    error(['G(3,3) indicates unrestrained heave. Use the submerged call ' ...
        'with T_1236 and delta_zeta_45.']);
end
if submergedCall && heaveIsRestored
    error(['G(3,3) indicates restored heave. Use the floating-vessel call ' ...
        'with kappa_126 and delta_zeta_345.']);
end

if ~isnumeric(aperiodicDamping) || ~isreal(aperiodicDamping) || ...
        any(~isfinite(aperiodicDamping(:)))
    error('The aperiodic damping parameters must be finite real numbers.');
end
if surfaceCall && any(aperiodicDamping(:) < 0)
    error('kappa_126 must contain three finite, nonnegative values.');
elseif submergedCall && any(aperiodicDamping(:) <= 0)
    error('T_1236 must contain four finite, positive time constants.');
end
if ~isnumeric(delta_zeta) || ~isreal(delta_zeta) || ...
        any(~isfinite(delta_zeta(:))) || any(delta_zeta(:) < 0)
    error('The damping-ratio increments must be finite and nonnegative.');
end
aperiodicDamping = reshape(aperiodicDamping, 1, []);
delta_zeta = reshape(delta_zeta, 1, []);

% Default plot flag
if nargin < 5 || isempty(plotFlag)
    plotFlag = 0;
end

% ------------------------------------------------------------------------------
%% Compute power-based equivalent matrices
% ------------------------------------------------------------------------------
freqs = vessel.freqs;

% Exclude artificial frequency omega = 10 rad/s representing infinity
idx = freqs < 10;
A = vessel.A(:,:,idx);
B = vessel.B(:,:,idx);
freqs = freqs(idx);

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

% PM wave spectrum parameters
constant = mssConstants();
alpha = 8.1e-3 * constant.g^2;
beta = 0.74;

% Initialize equivalent matrices
A_eq = zeros(6);
B_eq = zeros(6);

% Zero spectral moment m_0
S = alpha ./ freqs_fine.^5 .* exp(-beta * (omega_p ./ freqs_fine).^4);
m_0 = trapz(freqs_fine, S);

% Normalized wave spectrum
S_N = S / m_0;

% Loop over DOFs
for i = 1:6
    for j = 1:6
        A_ij = squeeze(A(i,j,:));
        B_ij = squeeze(B(i,j,:));

        A_interp = interp1(freqs, A_ij, freqs_fine, 'pchip');
        B_interp = interp1(freqs, B_ij, freqs_fine, 'pchip');

        A_eq(i,j) = trapz(freqs_fine, A_interp .* S_N);
        B_eq(i,j) = trapz(freqs_fine, B_interp .* S_N);
    end
end

% Power-based equivalent matrices
vessel.powerBased.omega_p = omega_p;
vessel.powerBased.A_eq = A_eq;
vessel.powerBased.B_eq = B_eq;

% ------------------------------------------------------------------------------
% System inertia matrix M = MRB + MA and restoring matrix G
% ------------------------------------------------------------------------------
vessel.MA = vessel.powerBased.A_eq; % Added mass matrix
vessel.M = vessel.MRB + vessel.MA;  % System inertia matrix
vessel.Minv = invQR(vessel.M);      % Inverse system inertia matrix
vessel.G = G;                       % Restoring matrix for heave, roll and pitch

% ------------------------------------------------------------------------------
% Viscous damping matrix Bv and total damping matrix D = B_eq + Bv
% ------------------------------------------------------------------------------
vessel.powerBased.Bv = zeros(6);
if surfaceCall
    % Floating vessel: relative viscous damping increments in DOFs 1, 2 and 6
    idx = [1 2 6];
    for k = 1:numel(idx)
        i = idx(k);
        vessel.powerBased.Bv(i,i) = aperiodicDamping(k) * ...
            vessel.powerBased.B_eq(i,i);
    end
    restoredIdx = [3 4 5];
else
    % Submerged vehicle: impose time constants in DOFs 1, 2, 3 and 6
    idx = [1 2 3 6];
    for k = 1:numel(idx)
        i = idx(k);
        targetDamping = vessel.M(i,i) / aperiodicDamping(k);
        viscousDamping = targetDamping - vessel.powerBased.B_eq(i,i);
        dampingTolerance = 1e-10 * max([1, abs(targetDamping), ...
            abs(vessel.powerBased.B_eq(i,i))]);
        if viscousDamping < -dampingTolerance
            error(['T_%d = %.4g s requires negative viscous damping. ' ...
                'Choose T_%d <= M(%d,%d)/B_eq(%d,%d) = %.4g s.'], ...
                i, aperiodicDamping(k), i, i, i, i, i, ...
                vessel.M(i,i) / vessel.powerBased.B_eq(i,i));
        end
        vessel.powerBased.Bv(i,i) = max(0, viscousDamping);
    end
    restoredIdx = [4 5];
end

% Restored modes: viscous damping-ratio increments
for k = 1:numel(restoredIdx)
    i = restoredIdx(k);
    vessel.powerBased.Bv(i,i) = 2 * delta_zeta(k) * ...
        sqrt(vessel.M(i,i) * vessel.G(i,i));
end

vessel.D = vessel.powerBased.B_eq + vessel.powerBased.Bv;

if surfaceCall
    vessel.powerBased.kappa_126 = aperiodicDamping;
    vessel.powerBased.delta_zeta_345 = delta_zeta;
else
    vessel.powerBased.T_1236 = aperiodicDamping;
    vessel.powerBased.delta_zeta_45 = delta_zeta;
end

vessel.powerBased.T_126 = zeros(1,3);
idx = [1 2 6];
for k = 1:3
    i = idx(k);
    vessel.powerBased.T_126(k) = vessel.M(i,i) / vessel.D(i,i);
end

% ------------------------------------------------------------------------------
%% Optional plotting
% ------------------------------------------------------------------------------
if plotFlag == 1
    plotAB_eq(vessel, 'A', 1);
    plotAB_eq(vessel, 'B', 1);
end

end
