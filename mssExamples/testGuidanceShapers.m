function testGuidanceShapers()
% testGuidanceShapers() checks the reference shapers of the GNC folder:
%
%   refModel         while the velocity limit holds v_d, the acceleration a_d
%                    must be the derivative of v_d (zero), not the unlimited
%                    value; this matters because a_d is the feed-forward term
%                    of the heading controllers (SIMotter, SIMmariner).
%   refModelPolyExp  a time vector is clipped element by element, t_max must be
%                    positive, and the header examples must run.
%   LOSobserver      K_f must be positive (K_f <= 0 gives a complex or frozen
%                    observer).
%
% The unsaturated path of refModel is compared with the third-order model of
% Fossen (2011), Eq. 10.27-10.28 (forward Euler), written out here.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eq. 10.27-10.28, p. 250 (third-order reference model)
%       and Eq. 12.138-12.142, p. 382 (feed-forward from the reference).
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
h = 0.01;

% ---- (1) refModel with the velocity limit active: step 30, v_max = 1
x = 0; v = 0; a = 0; n = round(300/h);
log = zeros(n, 3);
for k = 1:n
    log(k,:) = [x v a];
    [x, v, a] = refModel(x, v, a, 30, 1, 1.0, 0.5, h, 0);
end
onLimit = abs(log(2:end,2)) >= 1 - 1e-12;        % samples where v_d sits on v_max
dv = diff(log(:,2)) / h;                        % derivative of the stored velocity
errPlateau = max(abs(log(2:end,3) - dv));
nfail = nfail + report('a_d = dv_d/dt on the velocity limit, max error (rad/s^2)', errPlateau, 0.05);
nfail = nfail + report('time on the limit is not empty (s)', double(nnz(onLimit) < 100)*1, 0.5);
nfail = nfail + report('step reached, |x_d - 30| at 300 s', abs(log(end,1) - 30), 1e-3);
nfail = nfail + report('no overshoot, max(x_d) - 30', max(0, max(log(:,1)) - 30), 1e-2);

% ---- (2) refModel without saturation equals the third-order model (Euler)
x = 0; v = 0; a = 0; xr = 0; vr = 0; ar = 0; w = 1.0; z = 0.7; err = 0;
for k = 1:2000
    [x, v, a] = refModel(x, v, a, 1, 1e9, z, w, h, 0);
    jerk = -w^3*(xr - 1) - (2*z+1)*w^2*vr - (2*z+1)*w*ar;
    xr = xr + h*vr; vr = vr + h*ar; ar = ar + h*jerk;
    err = max(err, norm([x v a] - [xr vr ar]));
end
nfail = nfail + report('unsaturated refModel equals the third-order model, max error', err, 1e-12);

% ---- (3) the heading loop of SIMotter (case 1) on the matched plant
% tau_N = (T/K) a_d + (1/K) r_d - Kp (e + Td (r - r_d) + z/Ti), M r' + (M/T) r = tau_N
[~,~,M] = otter();
M66 = M(6,6); T = 1; K = T/M66; wn = 1.5; zeta = 1;
Kp = M66*wn^2; Kd = M66*(2*zeta*wn - 1/T); Td = Kd/Kp; Ti = 10/wn;
hc = 0.05; psi = 0; r = 0; zi = 0; pd = 0; rd = 0; ad = 0; emax = 0;
for k = 1:round(300/hc)
    ref = pi*(k*hc > 5);                         % 180 degree step at 5 s
    tau = (T/K)*ad + (1/K)*rd - Kp*(ssa(psi - pd) + Td*(r - rd) + zi/Ti);
    emax = max(emax, abs(ssa(psi - pd)));
    f = @(s) [s(2); (tau - (M66/T)*s(2))/M66];
    s0 = [psi; r];
    k1 = f(s0); k2 = f(s0 + hc/2*k1); k3 = f(s0 + hc/2*k2); k4 = f(s0 + hc*k3);
    s1 = s0 + hc/6*(k1 + 2*k2 + 2*k3 + k4);
    zi = zi + hc*ssa(psi - pd);
    psi = s1(1); r = s1(2);
    [pd, rd, ad] = refModel(pd, rd, ad, ref, deg2rad(10), 1.0, 1.0, hc, 1);
end
nfail = nfail + report('SIMotter heading loop, 180 deg step, max |psi - psi_d| (deg)', rad2deg(emax), 2.0);

% ---- (4) refModelPolyExp
t = [0 5 30];
y = refModelPolyExp(t, 20, 0, 10, 'quintic');
nfail = nfail + report('vector t, quintic, last element is x_final', abs(y(3) - 10), 1e-12);
yc = refModelPolyExp([-5 5 10 15 30], 20, 0, 10, 'cubic');
ycs = arrayfun(@(tt) refModelPolyExp(tt, 20, 0, 10, 'cubic'), [-5 5 10 15 30]);
nfail = nfail + report('vector t equals element-by-element scalar calls', max(abs(yc - ycs)), 1e-12);
ye = refModelPolyExp([-5 5], 20, 0, 10, 'exponential');
nfail = nfail + report('negative t is clipped to the start value', abs(ye(1) - 0), 1e-12);
nfail = nfail + check('t_max = 0 is refused', @() refModelPolyExp(5, 0, 0, 10, 'cubic'));
src = fileread(which('refModelPolyExp'));
bad = contains(src, '''qubic''');
nfail = nfail + report('header example names a method that exists', double(bad), 0.5);

% ---- (5) LOSobserver
nfail = nfail + check('K_f = 0 is refused', @() LOSobserver(0, 0, 1, 0.05, 0));
nfail = nfail + check('K_f < 0 is refused', @() LOSobserver(0, 0, 1, 0.05, -0.5));
[ang, rate] = LOSobserver(0, 0, 1, 0.05, 0.5);
nfail = nfail + report('K_f > 0 still works (finite, real)', double(~(isreal(ang) && isreal(rate) && isfinite(ang))), 0.5);

if nfail > 0
    error('testGuidanceShapers FAILED (%d check(s) wrong)', nfail);
end
disp('testGuidanceShapers PASSED');
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-68s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-68s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end

function bad = check(name, fcn)
try
    fcn(); fprintf('%-68s FAILED: no error\n', name); bad = 1;
catch
    fprintf('%-68s PASSED: error raised\n', name); bad = 0;
end
end
