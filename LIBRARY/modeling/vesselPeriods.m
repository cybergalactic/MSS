function [T,zeta,omega_d,omega_n] = vesselPeriods(w,MRB,Aw,Bw,C,method,verbose)
% [T,zeta,omega_d,omega_n] = vesselPeriods(w,MRB,Aw,Bw,C,method,verbose)
% computes the heave, roll, and pitch periods from frequency-dependent
% hydrodynamic data. Additionally, the undamped natural frequencies
% omega_n (rad/s), relative damping ratios zeta (-), and damped natural
% frequencies omega_d (rad/s) are computed from frequency-dependent added
% mass and damping data.
%
% Two formulations are supported:
%
%   1) Decoupled formulation ('decoupled'):
%      Each DOF is treated independently by solving the scalar implicit equation
%
%         f_i(omega_n) = omega_n^2 * ( MRB(ii) + Aii(omega_n) ) - C(ii) = 0
%
%      where Aii(omega) is obtained by interpolation of Aw(ii,ii,:) at omega.
%      The damping ratio is then computed as
%
%         zeta_i = Bii(omega_n) / ( 2*omega_n*( MRB(ii) + Aii(omega_n) ) )
%
%      and the damped natural frequency and period follow as
%
%         omega_d,i = omega_n,i * sqrt(1 - zeta_i^2)
%         T_i = 2*pi/omega_d,i
%
%   2) Coupled 6-DOF formulation ('coupled'):
%      The coupled undamped natural frequencies are obtained from the frequency-
%      dependent generalized eigenvalue problem
%
%         G x = lambda(omega) M(omega) x,   lambda = omega_n^2
%         M(omega) = MRB + A(omega),        G = C
%
%      where A(omega) is interpolated from Aw(:,:,k) at omega. For each selected
%      mode (heave/roll/pitch), the constraint is enforced using a scalar
%      root-finding method (fzero) on
%
%         g(omega) = omega^2 - lambda_sel(omega) = 0
%
%      where lambda_sel(omega) is the eigenvalue selected by modal participation
%      in the requested DOF (3,4,5).
%
%      Damping ratios are computed a posteriori at omega = omega_n using modal
%      (Rayleigh-quotient) quantities:
%
%         m_i = x_i' M(omega_n) x_i
%         k_i = x_i' G x_i
%         b_i = x_i' B(omega_n) x_i
%
%         zeta_i = b_i / (2*sqrt(k_i*m_i))
%         omega_d,i = omega_n,i * sqrt(1 - zeta_i^2)
%         T_i = 2*pi/omega_d,i
%
% INPUTS:
%   w       [N x 1]  frequency grid (rad/s), strictly increasing. The MSS
%                   infinity marker at 10 is excluded internally together
%                   with its Aw/Bw slices; a genuine zero sample is retained.
%   MRB     [6 x 6]  rigid-body mass/inertia matrix
%   Aw      [6 x 6 x N] added mass matrices on grid w
%   Bw      [6 x 6 x N] damping matrices on grid w (radiation + optional viscous)
%   C       [6 x 6]  restoring matrix (hydrostatic + optional Kp)
%   method  char     'decoupled' or 'coupled'
%   verbose logical  if true, print a small results table
%
% OUTPUTS (for DOFs 3,4,5 in this order):
%   T       [1 x 3]  periods (s) based on damped natural frequencies
%   zeta    [1 x 3]  relative damping ratios (-)
%   omega_d [1 x 3]  damped natural frequencies (rad/s)
%   omega_n [1 x 3]  undamped natural frequencies (rad/s)
%
% Reference:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., Wiley, Sec. 4.3.
%
% See also: fzero, eig, interp1, Hmtrx
%
% Author: Thor I. Fossen
% Date: 2026-02-13
% Revisions: 2026-10-27 Exclude the infinite-frequency marker from
%           interpolation and the natural-frequency search interval.
%   2026-10-29 Name the damped-frequency output and column omega_d.
%              Preserve the output order and undamped name omega_n.

% Natural periods require finite-frequency coefficients. In MSS, the
% radiation solution stored at w = 10 represents infinity, not A(10)/B(10).
w = w(:);
if size(Aw,3) ~= numel(w) || size(Bw,3) ~= numel(w)
    error('vesselPeriods:FrequencyDimensions', ...
        'Aw and Bw must have one matrix for each entry in w.');
end
finiteIdx = isfinite(w) & abs(w - 10) >= 10*eps(10);
w = w(finiteIdx);
Aw = Aw(:,:,finiteIdx);
Bw = Bw(:,:,finiteIdx);
if numel(w) < 2 || any(w < 0) || any(diff(w) <= 0)
    error('vesselPeriods:InvalidFrequencies', ...
        'At least two increasing, nonnegative finite frequencies are required.');
end

dofs  = [3 4 5];
w_min = w(1);
w_max = w(end);

% -------------------- Allocate outputs -----------------------------------
T       = zeros(1,3);
omega_d = zeros(1,3);
zeta    = zeros(1,3);
omega_n = zeros(1,3);

% -------------------- Main switch ----------------------------------------
switch lower(method)

    case 'decoupled'

        % --- Scalar fzero on each diagonal DOF (3,4,5) --------------------
        for k = 1:3
            i = dofs(k);

            % Interpolate added mass and damping (diagonal entries)
            Aii = @(w0) interp1(w, squeeze(Aw(i,i,:)), w0, 'pchip', 'extrap');
            Bii = @(w0) interp1(w, squeeze(Bw(i,i,:)), w0, 'pchip', 'extrap');

            % Undamped natural frequency: omega_n^2*(m + Aii(omega_n)) - Cii = 0
            f = @(wi) wi.^2 .* (MRB(i,i) + Aii(wi)) - C(i,i);
            omega_n(k) = fzero(f, [w_min w_max]);

            % Relative damping ratio: zeta = Bii(omega_n) / (2*omega_n*(m + Aii))
            meq = MRB(i,i) + Aii(omega_n(k));
            zeta(k) = Bii(omega_n(k)) / (2*omega_n(k)*meq);

            % Damped natural frequency and period
            omega_d(k) = omega_n(k) * sqrt(max(0, 1 - zeta(k)^2));
            T(k)     = 2*pi / omega_d(k);
        end

    case 'coupled'
        % --- Coupled 6-DOF: fzero on constraint omega^2 - lambda_sel(omega)=0
        for k = 1:3
            dof = dofs(k);

            % Root function: g(omega) = omega^2 - lambda_sel(omega)
            g = @(om) coupledConstraint(om, dof);

            % Solve for undamped natural frequency omega_n
            omega_n(k) = fzero(g, [w_min w_max]);

            % Modal vector and damping ratio at omega_n
            [x_mode, M_om, B_om] = coupledMode(omega_n(k), dof);

            m_modal = real(x_mode' * M_om * x_mode);
            k_modal = real(x_mode' * C    * x_mode);
            b_modal = real(x_mode' * B_om * x_mode);

            zeta(k)  = b_modal / (2*sqrt(k_modal*m_modal));
            omega_d(k) = omega_n(k) * sqrt(max(0, 1 - zeta(k)^2));
            T(k)     = 2*pi / omega_d(k);
        end

    otherwise
        error('Valid methods are ''decoupled'' or ''coupled''.');

end

% -------------------- Optional display -----------------------------------
if verbose
    fprintf('\n');
    fprintf('    DOF      T [s]      zeta     omega_d [rad/s]  omega_n [rad/s]\n');
    fprintf('    ---   ----------  --------  ---------------  ---------------\n');
    for k = 1:numel(dofs)
        fprintf('  %4d  %10.3f %10.3f   %10.3f     %10.3f\n', ...
            dofs(k), T(k), zeta(k), omega_d(k), omega_n(k));
    end
end


% ===== Nested helper functions =======================

function Aom = Aofw(om)
% Interpolate full added-mass matrix A(om) from tabulated Aw(:,:,k)
Aom = zeros(6,6);
for r = 1:6
    for c = 1:6
        Aom(r,c) = interp1(w, squeeze(Aw(r,c,:)), om, 'pchip', 'extrap');
    end
end
end

function Bom = Bofw(om)
% Interpolate full damping matrix B(om) from tabulated Bw(:,:,k)
Bom = zeros(6,6);
for r = 1:6
    for c = 1:6
        Bom(r,c) = interp1(w, squeeze(Bw(r,c,:)), om, 'pchip', 'extrap');
    end
end
end

function val = coupledConstraint(om, targetDof)
% Constraint function for coupled undamped natural frequency:
%   val(om) = om^2 - lambda_sel(om),
% where lambda_sel(om) is selected by modal participation in targetDof.
M = MRB + Aofw(om);
[V,D] = eig(C, M);
lam   = real(diag(D));
idx   = selectModeByParticipation(V, lam, targetDof);
val   = om^2 - lam(idx);
end

function [x_mode, M_om, B_om] = coupledMode(om, targetDof)
% Return selected eigenvector and interpolated matrices at frequency om.
M_om = MRB + Aofw(om);
B_om = Bofw(om);

[V,D] = eig(C, M_om);
lam   = real(diag(D));
idx   = selectModeByParticipation(V, lam, targetDof);

x_mode = V(:,idx);
% Normalize to avoid scaling issues (optional but harmless)
x_mode = x_mode / max(norm(x_mode), eps);
end

function idx = selectModeByParticipation(V, lam, dofSel)
% Select physically relevant mode by participation in a given DOF.
% Picks the eigenpair with positive eigenvalue and maximal |V(dofSel,i)|.
valid = find(isfinite(lam) & lam > 1e-9);
if isempty(valid)
    error('No positive finite eigenvalues found in coupled eigenproblem.');
end
[~,kk] = max(abs(V(dofSel,valid)));
idx    = valid(kk);
end

end % vesselPeriods
