function testNpsauvCrossFlowPitch()
% testNpsauvCrossFlowPitch() checks the cross-flow strip sum of npsauv.m.
%
% The strip at x = xL on the body x axis moves with the velocity
% nu_1 + omega x r = [u, v + r*xL, w - q*xL]', so the cross-flow load of the
% strip is F = -(rho/2) * drag_term * [0, vy, vz]'/Ucf * dxL and its moment
% about the origin is r x F = [0, -xL*Fz, xL*Fy]'. npsauv.m already uses
% w_r - xL*q in Cz and in drag_term; the pitch moment M = -xL*Fz then needs
% (xL*q - w_r) in Cm, since tau_crossflow(5) = -(rho/2)*Cm. The file had
% (w_r + xL*q), which agrees with the physical load only when w_r = 0.
%
% The test takes the strip loop of npsauv.m as it is written in the file (the
% lines from "dxL = L/10" to "tau_crossflow = ..."), runs it on 203 relative
% velocities (the three of the review, 200 seeded ones) and compares the
% vector tau_crossflow with the strip sum built above from the rigid-body
% point velocity and r x F. A sign error in the pitch term fails the pitch
% moment, the other five components pass on both versions. The agreement of
% the fixed file is 5e-7 of the largest load, not machine precision, because
% npsauv.m adds epsilon_drag = 1e-6 to Ucf.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eq. 2.14 (velocity of a point of a rigid body) and
%       Section 6.6 (strip theory).
%
% Author:    Enio Krizman
% Date:      2 November 2026

% Constants as in npsauv.m (the loop reads them from the workspace)
L = 5.3; Cdy = 0.5; Cdz = 0.6; rho = mssConstants().rho_water; %#ok<NASGU>

% The strip loop of npsauv.m, copied from the file
txt = fileread(which('npsauv'));
blk = regexp(txt, 'dxL\s*=\s*L/10.*?tau_crossflow\s*=[^\n]*', 'match', 'once');
if isempty(blk)
    error('testNpsauvCrossFlowPitch: strip loop not found in npsauv.m');
end
v_r = 0; w_r = 0; q = 0; r = 0;
eval(blk);   % defines dxL, Hx, Bx and tau_crossflow from v_r, w_r, q, r

rng(1);
cases = [0.2 0.1; 0.5 0.2; 0.3 -0.1; 0.3*randn(200, 2)]';   % columns [w_r; q]
errAll = 0; errPitch = 0; scale = 0;
for k = 1:size(cases, 2)
    v_r = 0; r = 0;
    if k > 3, v_r = 0.3 * randn; r = 0.1 * randn; end
    w_r = cases(1, k); q = cases(2, k);
    eval(blk);
    tau_file = tau_crossflow;

    tau_ref = zeros(6, 1);
    for xL = -L/2:dxL:L/2
        vy = v_r + r * xL;  vz = w_r - q * xL;       % point velocity
        Ucf = sqrt(vy^2 + vz^2);
        if Ucf == 0, continue; end
        drag = Cdy * Hx * vy^2 + Cdz * Bx * vz^2;
        F = -(rho/2) * drag * [0; vy; vz] / Ucf * dxL;
        tau_ref = tau_ref + [F; cross([xL; 0; 0], F)];
    end
    e = abs(tau_file - tau_ref);
    errAll = max(errAll, max(e));
    errPitch = max(errPitch, e(5));
    scale = max(scale, max(abs(tau_ref)));
end
if k > 0 && ~(errPitch <= 1e-6 * scale)
    fprintf('pitch moment: max |file - strip sum| = %.4g N m (largest load %.4g)\n', errPitch, scale);
end
fprintf('three cases of the review, pitch moment in the file (N m): ');
for k = 1:3
    w_r = cases(1, k); q = cases(2, k); v_r = 0; r = 0; eval(blk);
    fprintf('%.4g ', tau_crossflow(5));
end
fprintf('\n');

if errAll > 1e-6 * scale
    error('testNpsauvCrossFlowPitch FAILED (max |tau_crossflow - strip sum| = %.4g, pitch %.4g)', ...
        errAll, errPitch);
end
fprintf('npsauv cross-flow strip sum vs rigid-body strip sum: max error %.2e\n', errAll);
disp('testNpsauvCrossFlowPitch PASSED');
end
