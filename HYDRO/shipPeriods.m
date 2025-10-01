function [w_345,T_345] = shipPeriods(vessel, dof)
% shipPeriods  Natural frequencies (rad/s) and periods (s) in heave/roll/pitch
%
%   [w_345,T_345] = shipPeriods(vessel, dof)
%
% Computes the natural frequencies (rad/s) of a 6-DOF vessel model using
% both coupled and decoupled formulations:
%
%   - Coupled case: The problem is solved in the CG by iteratively
%     evaluating the generalized eigenvalue problem
%
%           | G − ω² ( M_RB + A(ω) ) | = 0
%
%     The eigenvalues ω² yield the coupled natural frequencies.
%
%   - Decoupled case: The rigid-body, added mass, and restoring matrices
%     are transformed from the CO to the CG using Hmtrx.m. The natural 
%     requencies in heave,roll, and pitch are then found by solving the scalar 
%     implicit equation
%
%           ω_i = sqrt( G_CG(i,i) / ( M_RB_CG(i,i) + A_CG(i,i) ) )
%
% Natural periods follow as  T_i = 2π / ω_i.
%
% Inputs (vessel struct expressed in the CO (centerline in the WL, midships):
%   vessel.MRB   : [6x6]   rigid-body mass/inertia at BODY/CO
%   vessel.A     : [6x6xN] added mass at BODY/CO on grid vessel.freqs
%   vessel.B     : [6x6xN] potential damping on same grid (not used here)
%   vessel.C     : [6x6xN] hydrostatic restoring at BODY/CO (use C(:,:,1))
%   vessel.freqs : [N x 1] monotonically increasing frequency grid [rad/s]
%   vessel.CG    : [1x3]   position of CG expressed in BODY (r_CO→CG)
%
% Mode selector:
%   dof = 0 : Decoupled (CG), compute heave/roll/pitch 
%   dof = 1 : Coupled (CO), compute heave/roll/pitch 
%
% Outputs:
%   w_345 (3x1) : Natural frequencies [ w3; w4; w5 ] in rad/s
%   T_345 (3x1) : Natural periods     [ T3; T4; T5 ] in s
%
% Notes:
%   • Coupled solution uses matrices at the BODY/CO as given.
%   • Decoupled solution transforms MRB, A(ω), C from CO→CG using Hmtrx(r),
%     then solves ω = sqrt( C_ii / ( MRB_ii + A_ii(ω) ) ) for i∈{3,4,5}.
%   • vessel.B is accepted for completeness but not used for undamped ω_n.

% Vessel matrices expressed in the CO
MRB   = vessel.MRB;
A     = vessel.A;
G     = vessel.C(:,:,1);
freqs = vessel.freqs(:);

% [x y z] of CG in BODY (from CO)
r_CG = [ vessel.main.CG(1)
         0
         vessel.main.T - vessel.main.CG(3) ];

N = numel(freqs);

% Screw transformation from the CO to the CG
H    = Hmtrx(r_CG);
Hinv = inv(H);

MRB_CG = Hinv' * MRB  * Hinv;
G_CG   = Hinv' * G * Hinv;

% Transform A(ω) to CG for all grid points
A_CG = zeros(6,6,N);
for k = 1:N
    A_CG(:,:,k) = Hinv' * A(:,:,k) * Hinv;
end

% decoupled scalar solver (in CG)
decoupled_one = @(idxCG) decoupled_scalar(idxCG, MRB_CG, G_CG, A_CG, freqs);

% indices for heave/roll/pitch
idxMap = [3 4 5];

if dof == 1
    % Coupled 6-DOF 
    % 1) Build decoupled guesses in CG (robust initial guesses)
    w_dec = zeros(3,1);
    
    for k = 1:3
        w_dec(k) = decoupled_one(idxMap(k));
    end

    % 2) Iterate generalized eigenproblem in CO using those guesses
    w_cpl = zeros(3,1);
    w_cpl(1) = mode_iter_weighted(3, w_dec(1), MRB, G, A, freqs); % heave
    w_cpl(2) = mode_iter_weighted(4, w_dec(2), MRB, G, A, freqs); % roll
    w_cpl(3) = mode_iter_weighted(5, w_dec(3), MRB, G, A, freqs); % pitch
    w_345 = w_cpl(:);

else
    % Decoupled in Cthe G 
    w_dec = zeros(3,1);
    for k = 1:3
        w_dec(k) = decoupled_one(idxMap(k));
    end
    w_345 = w_dec(:);

end

% Periods
T_345 = 2*pi ./ w_345;

end

% ==============================================================================
% Scalar decoupled fixed-point for a single oscillatory DOF (in CG)
% ==============================================================================
function w = decoupled_scalar(idxCG, MRB_CG, G_CG, A_CG, freqs)
k    = G_CG(idxCG, idxCG);
Aii  = squeeze(A_CG(idxCG, idxCG, :));

% Initial guess from mid-band diagonal
Amid = A_CG(:,:,max(1,round(numel(freqs)/2)));
w0   = sqrt( max(k,0) / max(MRB_CG(idxCG,idxCG) + Amid(idxCG,idxCG), eps) );

w     = max(w0, 1e-6);
tol   = 1e-10;
maxit = 60;

for it = 1:maxit
    wc    = min(max(w, freqs(1)), freqs(end));
    Aii_w = interp1(freqs, Aii, wc, 'pchip');
    w_new = sqrt( k / max(MRB_CG(idxCG,idxCG) + Aii_w, eps) );

    if ~isfinite(w_new) || w_new <= 0, break; end
    
    if abs(w_new - w) / max(w,1e-6) < tol
        w = w_new; break;
    end
    
    w = w_new;
end

end

% ==============================================================================
% Interpolate full A(ω) at a scalar frequency
% ==============================================================================
function Aint = Aofw(w, Acoeff, freqs)
wc = min(max(w, freqs(1)), freqs(end));
Aint = zeros(6,6);

for r = 1:6
    for c = 1:6
        Aint(r,c) = interp1(freqs, squeeze(Acoeff(r,c,:)), wc, 'pchip');
    end
end

end

% ==============================================================================
% Coupled mode finder with participation + proximity weighting (in CO)
% ==============================================================================
function wstar = mode_iter_weighted(target, wguess, MRB_CO, G_CO, Acoeff_CO, freqs)
w     = max(wguess, 1e-6);
tol   = 1e-6;
maxit = 60;
wstar = wguess;                         % Fallback
sigma = max(0.35*wguess, 0.2);          % Proximity bandwidth [rad/s]

for it = 1:maxit
    M = MRB_CO + Aofw(w, Acoeff_CO, freqs);
    [V,D] = eig(G_CO, M);
    lam   = real(diag(D));
    pos   = find(isfinite(lam) & lam > 1e-9);
    if isempty(pos), return; end

    wcand = sqrt(lam(pos));   wcand = wcand(:);
    P     = abs(V(target, pos)); P = P(:);
    if isempty(wcand) || all(~isfinite(wcand)) || all(P==0), return; end

    S = (P./max(P)) .* exp(-((wcand - wguess).^2)/(sigma^2));
    S(~isfinite(S)) = -Inf;

    [~,idx] = max(S);
    w_new   = wcand(idx);

    if ~isfinite(w_new) || w_new <= 0, return; end
    
    if abs(w_new - w) / max(w,1e-6) < tol
        wstar = w_new; return;
    end
    
    w = w_new; wstar = w;
end

end