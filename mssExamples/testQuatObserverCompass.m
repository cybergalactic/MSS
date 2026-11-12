function testQuatObserverCompass()
% testQuatObserverCompass() checks the compass form of quatObserver.m
% (imu_meas = [f' w' psi]). With exact measurements and a constant attitude
% the observer must stay at the true attitude. The injection compared the BODY
% x-axis with the NED heading vector rotated into the BODY frame; the two agree
% only at zero roll and pitch, so at a constant tilt the observer converged to
% an attitude that was wrong by about k2*sin(theta)/(k1+k2) in pitch (0.32 deg
% at 10 deg pitch, 1.32 deg at 45 deg for the gains below). The compass gives a
% rotation about the NED down axis, which the injection must use.
%
% Also checked: the observer still converges from an initial attitude error
% (level and tilted), and the magnetometer form is unchanged.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   R. Mahony, T. Hamel, J.-M. Pflimlin (2008). Nonlinear complementary
%       filters on the special orthogonal group. IEEE TAC 53(5), Eq. 32a-c.
%   H. F. Grip, T. I. Fossen, T. A. Johansen, A. Saberi (2011). Nonlinear
%       observer for GNSS-aided inertial navigation with quaternion-based
%       attitude estimation, Eq. 2-3.
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
g_n = [0 0 9.81]'; k1 = 1500; k2 = 50; Ki = 0.1*eye(3); h = 0.001;

% ---- constant tilt, exact measurements, started at the true attitude
for att = [0 10 0 40; 0 20 0 40; 0 45 0 40; 20 20 0 40; -30 15 0 -130]'
    phi = deg2rad(att(1)); th = deg2rad(att(2)); psi = deg2rad(att(4));
    qt = euler2q(phi, th, psi); q = qt; b = zeros(3,1);
    R = Rzyx(phi, th, psi); f = (-R' * g_n)';
    for i = 1:8000
        [q, b] = quatObserver(q, b, h, Ki, k1, k2, [1 0 0]', [f 0 0 0 psi]);
    end
    nfail = nfail + report(sprintf('roll %3d, pitch %2d deg: attitude error (deg)', att(1), att(2)), ...
        attErr(q, qt), 1e-3);
end

% ---- convergence from an initial error, tilted attitude, with a gyro bias
phi = deg2rad(10); th = deg2rad(30); psi = deg2rad(40);
qt = euler2q(phi, th, psi); R = Rzyx(phi, th, psi); f = (-R' * g_n)';
e0 = deg2rad([5 -4 8])'; q = quatprod(qt, expLocal(-e0)); b = zeros(3,1);
bg = [0.01 -0.02 0.015];
for i = 1:20000
    [q, b] = quatObserver(q, b, h, Ki, k1, k2, [1 0 0]', [f bg psi]);
end
nfail = nfail + report('from a 10 deg initial error with a gyro bias: attitude error (deg)', attErr(q, qt), 1e-2);
fprintf('  (gyro bias estimate error after 20 s: %.2e rad/s of %.2e)\n', norm(b - bg'), norm(bg));

% ---- magnetometer form unchanged (exact at a tilted attitude)
mref = [20000; 1500; 45000];
qt = euler2q(deg2rad(20), deg2rad(30), deg2rad(40)); R = Rquat(qt); f = (-R' * g_n)'; m = (R' * mref)';
q = qt; b = zeros(3,1);
for i = 1:4000
    [q, b] = quatObserver(q, b, h, Ki, 1000, 500, mref, [f 0 0 0 m]);
end
nfail = nfail + report('magnetometer form at a tilted attitude: attitude error (deg)', attErr(q, qt), 1e-3);

if nfail > 0
    error('testQuatObserverCompass FAILED (%d check(s) wrong)', nfail);
end
disp('testQuatObserverCompass PASSED');
end

function e = attErr(q, qt)
dq = quatprod([qt(1); -qt(2:4)], q);
e = rad2deg(2*norm(dq(2:4)));
end

function q = expLocal(th)
a = norm(th);
if a < 1e-12, q = [1; th/2]; else, q = [cos(a/2); sin(a/2)*th/a]; end
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-70s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-70s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
