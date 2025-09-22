function [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G)
% [T1,T2,T6,w3,w4,w5,zeta3,zeta4,zeta5] = modalFromMDG(M,D,G)
%
% Compute modal properties of a 6-DOF vessel model:
%   M nu_dot + D nu+ G eta = 0
%
% Uses eigenvalues/eigenvectors of the state matrix to map:
%   - Real poles  -> surge (1), sway (2), yaw (6)
%   - Complex pairs -> heave (3), roll (4), pitch (5)

    n = size(M,1);
    if n ~= 6, error('Expected 6x6 matrices'); end

    % Build state matrix: x=[q; qd], xdot = A x
    A = [zeros(n) eye(n); -M\G  -M\D];

    [V, Lambda] = eig(A);
    ev = diag(Lambda);

    % ----- Real (aperiodic) modes for surge,sway,yaw -----
    idx_real = find(abs(imag(ev)) < 1e-8 & real(ev) < 0);
    p_real   = real(ev(idx_real));
    Vq_real  = V(1:n, idx_real);

    dofs_aps = [1 2 6]; % surge,sway,yaw
    T = nan(1,3);
    magMat = abs(Vq_real(dofs_aps,:));
    for pass = 1:3
        [iRow,iCol] = find(magMat == max(magMat(:)),1);
        T(iRow) = -1/p_real(iCol);
        magMat(iRow,:) = -inf;
        magMat(:,iCol) = -inf;
    end
    T1 = T(1); T2 = T(2); T6 = T(3);

    % ----- Complex (oscillatory) modes for heave,roll,pitch -----
    idx_cplx = find(imag(ev) > 0);
    ev_c     = ev(idx_cplx);
    Vq_c     = V(1:n, idx_cplx);

    sigma = real(ev_c);
    omegad = imag(ev_c);
    wn_modes   = sqrt(sigma.^2 + omegad.^2);
    zeta_modes = -sigma ./ wn_modes;

    dofs_osc = [3 4 5]; % heave, roll, pitch
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