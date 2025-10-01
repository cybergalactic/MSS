function [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G)
% [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G)
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
%   • Real negative poles  → surge (DOF 1), sway (DOF 2), yaw (DOF 6)
%       - Identified as aperiodic modes with decay constants.
%       - Time constants Ti = -1/p_real.
%
%   • Complex conjugate poles → heave (DOF 3), roll (DOF 4), pitch (DOF 5)
%       - Identified as oscillatory modes.
%       - Natural frequencies:  wi = sqrt(σ² + ωd²)   [rad/s]
%       - Damping ratios:       ζi = -σ / wi          [-]
%
% Outputs:
%   T1,T2,T6         : Time constants in surge, sway, and yaw [s]
%   w3,w4,w5         : Natural frequencies in heave, roll, pitch [rad/s]
%   zeta3,zeta4,zet5 : Relative damping ratios in heave, roll, pitch [-]
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

n = size(M,1);
if n ~= 6, error('Expected 6x6 matrices'); end

% Build state matrix: x=[eta; nu], xdot = A x
A = [zeros(n) eye(n); -M\G  -M\D];

[V, Lambda] = eig(A);
ev = diag(Lambda);

% Real (aperiodic) modes for surge,sway,yaw
idx_real = find(abs(imag(ev)) < 1e-8 & real(ev) < 0);
p_real   = real(ev(idx_real));
Vq_real  = V(1:n, idx_real);

dofs_aps = [1 2 6]; % Surge,sway,yaw
T = nan(1,3);
magMat = abs(Vq_real(dofs_aps,:));
for pass = 1:3
    [iRow,iCol] = find(magMat == max(magMat(:)),1);
    T(iRow) = -1/p_real(iCol);
    magMat(iRow,:) = -inf;
    magMat(:,iCol) = -inf;
end
T1 = T(1); T2 = T(2); T6 = T(3);

% Complex (oscillatory) modes for heave,roll,pitch 
idx_cplx = find(imag(ev) > 0);
ev_c     = ev(idx_cplx);
Vq_c     = V(1:n, idx_cplx);

sigma = real(ev_c);
omegad = imag(ev_c);
wn_modes   = sqrt(sigma.^2 + omegad.^2);
zeta_modes = -sigma ./ wn_modes;

dofs_osc = [3 4 5]; % Heave, roll, pitch
w = nan(1,3);
z = nan(1,3);
magMat = abs(Vq_c(dofs_osc,:));
for pass = 1:3
    [iRow,iCol] = find(magMat == max(magMat(:)),1);
    w(iRow) = wn_modes(iCol);
    z(iRow) = zeta_modes(iCol);
    magMat(iRow,:) = -inf;
    magMat(:,iCol) = -inf;
end
w3 = w(1); w4 = w(2); w5 = w(3);
zeta3 = z(1); zeta4 = z(2); zeta5 = z(3);
end