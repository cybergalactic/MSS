function testRelativeVelocityTerms()
% testRelativeVelocityTerms() checks that npsauv.m and osv.m use the velocity
% relative to the water, nu_r = nu - nu_c, in the two terms that used the
% velocity over ground.
%
% In a uniform current that is constant in NED every force on a craft
% depends on the motion relative to the water (Fossen 2011, Eqs. 8.142-8.147).
% For a state x_r = [nu_r; eta] in still water and the same state moving with
% the current, x = [nu_r + nu_c; eta], the accelerations must satisfy
%
%   xdot(1:6) - xdot_r(1:6) = nu_c_dot,
%
% with nu_c and nu_c_dot from currentVelocity.m (the full-attitude form).
% testCurrentKinematics.m tests this where the two terms have no effect
% (npsauv with the propeller stopped, osv drifting with the current). Here
% they are active:
%   npsauv: the propeller advance ratio 0.012*n_p/u and the slip epsilon,
%           with the propeller running;
%   osv:    the added-mass Coriolis matrix, CA = m2c(MA, nu) instead of
%           m2c(MA, nu_r) (Fossen 2011, Eq. 8.136), with the vessel moving
%           through the water.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Section 8.3, pp. 221-225, Eq. 8.136.
%
% Author:    Enio Krizman
% Date:      2 November 2026

N = 200;
tol = 1e-10;
rng(2);
names = {'npsauv (propeller running)', 'osv (moving through the water)'};
failed = false;

for i = 1:numel(names)
    maxError = 0;
    for k = 1:N
        Vc = rand; betaVc = 2*pi*(rand - 0.5);
        Theta = [deg2rad(15)*(2*rand(2,1) - 1); pi*(2*rand - 1)];
        pos = [100*randn(2,1); 20*rand];
        nu_r = [1.5*rand + 0.2; 0.3*randn(2,1); 0.3*randn(3,1)];
        if i == 1
            w_c = 0.2*(2*rand - 1);
            delta = deg2rad(20)*(2*rand(4,1) - 1);
            n_p = 200 + 1200*rand;                    % running propeller
            x_r = [nu_r; pos; Theta; delta; n_p];
            ui = [deg2rad(20)*(2*rand(4,1) - 1); n_p];
            still = @(x) npsauv(x, ui);
            moving = @(x) npsauv(x, ui, Vc, betaVc, w_c);
            att = 10:12;
        else
            w_c = 0;
            pos(3) = 0.1*randn;
            x_r = [nu_r; pos; Theta];
            ui = [2*rand(4,1) - 1; pi*(2*rand(2,1) - 1)];
            still = @(x) osv(x, ui, 0, 0);
            moving = @(x) osv(x, ui, Vc, betaVc);
            att = 10:12;
        end
        [nu_c, nu_c_dot] = currentVelocity(Vc, betaVc, w_c, Theta, nu_r(4:6));
        x = x_r; x(1:6) = x_r(1:6) + nu_c;
        xdot_r = still(x_r);
        xdot = moving(x);
        e = max(abs(xdot(1:6) - xdot_r(1:6) - nu_c_dot)) ...
            / max(1, max(abs(xdot_r(1:6))));
        maxError = max(maxError, e);
        maxError = max(maxError, max(abs(xdot(att) - xdot_r(att))));
    end
    if maxError < tol, status = 'PASSED'; else, status = 'FAILED'; failed = true; end
    fprintf('%-34s %s: max error %.2e\n', names{i}, status, maxError);
end

if failed
    error('testRelativeVelocityTerms FAILED (tolerance %.0e)', tol);
end
fprintf('testRelativeVelocityTerms PASSED (tolerance %.0e)\n', tol);
end
