function testCrossFlowDragStripWidth()
% testCrossFlowDragStripWidth() checks the strip widths of crossFlowDrag.m.
%
% A strip of the hull at [x,0,0] moves with the velocity
% nu_1 + omega x r = [u, v + r x, w - q x]'. The transverse flow meets the
% strip with the width projected normal to it: the draft T for the sway
% flow, the beam B for the heave flow. With the 2-D coefficient Cd of the
% selected model the load per strip is
%
%    F = -1/2 rho Cd dx [0, T |vy| vy, B |vz| vz]'    and    M = r x F.
%
% The file used T for both, which is right only for B = T (a circular
% cylinder, the AUV case). The test builds the 20-strip sum above, from the
% rigid-body point velocity and r x F, and compares the full vector
% tau_crossflow for three surface hull sections (B > T, B < T, B = T) and
% for the cylinder model of the AUV (B = T = D), 200 seeded relative
% velocities each. The sway and yaw components agree on both versions, the
% heave force and pitch moment fail for B ~= T before the fix.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eq. 2.14 and Section 6.6 (strip theory).
%
% Author:    Enio Krizman
% Date:      2 November 2026

rho = mssConstants().rho_water;
hulls = { 'Hoerner',  60, 12.0, 3.0      % B > T, a ship hull
          'Hoerner',  60,  3.0, 6.0      % B < T
          'Hoerner',  40,  4.0, 4.0      % B = T
          'cylinder', 1.6, 0.19, 0.19 }; % AUV: B = T = D
rng(1);
fail = 0;
for h = 1:size(hulls, 1)
    model = hulls{h,1}; L = hulls{h,2}; B = hulls{h,3}; T = hulls{h,4};
    errAll = zeros(6,1); scale = zeros(6,1);
    for k = 1:200
        nu_r = [1; 0.3*randn(2,1); 0.1*randn(3,1)];
        nu_r(1) = 1 + randn;
        tau_file = crossFlowDrag(L, B, T, nu_r, model);

        % Reference: 20 midpoint strips, point velocity, r x F
        switch model
            case 'Hoerner',  Cd = Hoerner(B, T);
            case 'cylinder', Cd = cylinderDrag(L, B, nu_r);
        end
        dx = L/20; tau_ref = zeros(6,1);
        for i = 1:20
            x  = -L/2 + (i - 0.5)*dx;
            vp = nu_r(1:3) + cross(nu_r(4:6), [x; 0; 0]);   % point velocity
            F  = -0.5*rho*Cd*dx * [0; T*abs(vp(2))*vp(2); B*abs(vp(3))*vp(3)];
            tau_ref = tau_ref + [F; cross([x; 0; 0], F)];
        end
        errAll = max(errAll, abs(tau_file - tau_ref));
        scale  = max(scale, abs(tau_ref));
    end
    rel = max(errAll ./ max(scale, eps));
    ok = rel < 1e-12;
    fprintf('%-8s L=%5.1f B=%5.2f T=%5.2f : max relative error %9.2e  (heave force %9.2e, pitch moment %9.2e)  %s\n', ...
        model, L, B, T, rel, errAll(3)/scale(3), errAll(5)/scale(5), passFail(ok));
    fail = fail + ~ok;
end

if fail > 0
    error('testCrossFlowDragStripWidth: %d case(s) FAILED', fail);
end
disp('testCrossFlowDragStripWidth PASSED');
end

function s = passFail(ok)
if ok, s = 'ok'; else, s = 'FAIL'; end
end
