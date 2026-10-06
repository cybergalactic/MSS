function [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5,T3] = modalFromMDG(M,D,G,printFlag)
% [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G)
% [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G,printFlag)
% [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5,T3] = modalFromMDG(M,D,G,printFlag)
%
% Compute modal properties of a linear 6-DOF vessel model:
%
%     eta_dot = nu
%     M * nu_dot + D * nu + G * eta = 0
%
% where
%   M : 6x6 inertia + added mass matrix
%   D : 6x6 linear damping matrix
%   G : 6x6 hydrostatic restoring stiffness matrix
%
% The method forms the state-space system
%
%     x = [eta; nu],   xdot = A * x
%     A = [  0    I ;
%           -M\G -M\D ]
%
% and computes the eigenvalues/vectors of A.
%
% Mapping of eigenmodes:
%   • Real negative poles  → surge (DOF 1), sway (DOF 2), yaw (DOF 6),
%       and heave (DOF 3) for a submerged vehicle with G(3,3) = 0
%       - Identified as aperiodic modes with decay constants.
%       - Time constants Ti = -1/p_real.
%
%   • Complex conjugate poles → heave (DOF 3) for a floating vessel,
%       roll (DOF 4), and pitch (DOF 5)
%       - Identified as oscillatory modes.
%       - Natural frequencies:  wi = sqrt(σ² + ωd²)   [rad/s]
%       - Damping ratios:       ζi = -σ / wi          [-]
%
% Outputs:
%   T1,T2,T6         : Time constants in surge, sway, and yaw [s]
%   w3,w4,w5         : Natural frequencies in heave, roll, pitch [rad/s].
%                      For a submerged vehicle, w3 is NaN.
%   zeta3,zeta4,zeta5: Relative damping ratios in heave, roll, pitch [-].
%                      For a submerged vehicle, zeta3 is NaN.
%   T3               : Heave time constant [s] for a submerged vehicle.
%                      For a floating vessel, T3 is NaN (optional output).
%
% Optional input:
%   printFlag         : Print the modal properties when true (default false)
%                       Example: modalFromMDG(M,D,G,1)
%
% Notes:
%   - Assumes strictly 6 DOFs: surge (1), sway (2), heave (3),
%     roll (4), pitch (5), yaw (6).
%   - Largest modal participation is used to assign each mode
%     to its corresponding DOF.
%   - Purely real modes correspond to low-frequency rigid-body
%     motions; oscillatory modes correspond to hydrostatic restoring.
%
% Author:    Thor I. Fossen
% Date:      2025-10-01
% Revisions:
%   2026-10-04 Added optional printing of modal properties.
%              Added submerged-vehicle modes using the G(3,3) check.

if nargin < 4
    printFlag = false;
end
if ~(isscalar(printFlag) && (islogical(printFlag) || isnumeric(printFlag)))
    error('printFlag must be a logical or numeric scalar');
end

n = size(M,1);
if n ~= 6, error('Expected 6x6 matrices'); end

% Build state matrix: x=[eta; nu], xdot = A x
A = [zeros(n) eye(n); -M\G  -M\D];

[V, Lambda] = eig(A);
ev = diag(Lambda);

% G(3,3) distinguishes restored heave for floating vessels from aperiodic
% heave for fully submerged vehicles.
tolG = 1e-10 * max(1, norm(G, 'fro'));
heaveIsRestored = abs(G(3,3)) > tolG;

% Real (aperiodic) modes
idx_real = find(abs(imag(ev)) < 1e-8 & real(ev) < 0);
p_real   = real(ev(idx_real));
Vq_real  = V(1:n, idx_real);

if heaveIsRestored
    dofs_aps = [1 2 6];   % Floating vessel: surge, sway, yaw
else
    dofs_aps = [1 2 3 6]; % Submerged vehicle: surge, sway, heave, yaw
end
T = nan(1,6);
magMat = abs(Vq_real(dofs_aps,:));
for pass = 1:min(size(magMat))
    [maxParticipation, linearIndex] = max(magMat(:));
    if isempty(maxParticipation) || ~isfinite(maxParticipation)
        break
    end
    [iRow,iCol] = ind2sub(size(magMat), linearIndex);
    T(dofs_aps(iRow)) = -1/p_real(iCol);
    magMat(iRow,:) = -inf;
    magMat(:,iCol) = -inf;
end
T1 = T(1); T2 = T(2); T3 = T(3); T6 = T(6);

% Complex (oscillatory) modes for heave, roll, pitch
idx_cplx = find(imag(ev) > 0);
ev_c     = ev(idx_cplx);
Vq_c     = V(1:n, idx_cplx);

sigma = real(ev_c);
omegad = imag(ev_c);
wn_modes   = sqrt(sigma.^2 + omegad.^2);
zeta_modes = -sigma ./ wn_modes;

if heaveIsRestored
    dofs_osc = [3 4 5]; % Floating vessel: heave, roll, pitch
else
    dofs_osc = [4 5];   % Submerged vehicle: roll, pitch
end
w = nan(1,6);
z = nan(1,6);
magMat = abs(Vq_c(dofs_osc,:));
for pass = 1:min(size(magMat))
    [maxParticipation, linearIndex] = max(magMat(:));
    if isempty(maxParticipation) || ~isfinite(maxParticipation)
        break
    end
    [iRow,iCol] = ind2sub(size(magMat), linearIndex);
    w(dofs_osc(iRow)) = wn_modes(iCol);
    z(dofs_osc(iRow)) = zeta_modes(iCol);
    magMat(iRow,:) = -inf;
    magMat(:,iCol) = -inf;
end
w3 = w(3); w4 = w(4); w5 = w(5);
zeta3 = z(3); zeta4 = z(4); zeta5 = z(5);

if printFlag
    fprintf('%s\n','');
    fprintf('%s\n','------------------------------------------------------------------');
    fprintf('MODAL PROPERTIES OF THE VESSEL MODEL\n');
    fprintf('%s\n','------------------------------------------------------------------');
    fprintf('Aperiodic mode      Time constant [s]\n');
    fprintf('Surge, T1                 %10.4f\n', T1);
    fprintf('Sway,  T2                 %10.4f\n', T2);
    if ~heaveIsRestored
        fprintf('Heave, T3                 %10.4f\n', T3);
    end
    fprintf('Yaw,   T6                 %10.4f\n', T6);
    fprintf('\nOscillatory mode    Frequency [rad/s]    Damping ratio [-]\n');
    if heaveIsRestored
        fprintf('Heave, w3                 %10.4f            %10.4f\n', w3, zeta3);
    end
    fprintf('Roll,  w4                 %10.4f            %10.4f\n', w4, zeta4);
    fprintf('Pitch, w5                 %10.4f            %10.4f\n\n', w5, zeta5);
end
end
