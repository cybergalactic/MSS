function testInsBiasDecay()
% testInsBiasDecay() checks that the INS Kalman functions decay the nominal
% sensor biases with the time constant of their error-state model. The error
% model of ins_mekf, ins_mekf_psi, ins_euler, ins_ahrs, ins_heave and quatMEKF
% contains the Gauss-Markov term  d(delta_b)/dt = -delta_b/T  but the nominal
% bias b was held between corrections, so the estimate b + delta_b does not
% follow the process that the filter assumes (db/dt = -b/T + w). The decay is
% exact for that model:
%
%   b[k+1] = exp(-h/T) b[k]
%
% Without an aiding measurement (and with a diagonal covariance) the biases of
% the output state must therefore equal exp(-h/T) times those of the input.
%
% The last check runs the filter for 600 s without aiding: the gyro bias
% follows exp(-t/T) b0.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   Sola, J. (2017). Quaternion kinematics for the error-state Kalman filter.
%       Eq. 235 (bias error dynamics, the random-walk case of the Gauss-Markov
%       model used here).
%   Bar-Shalom, Li, Kirubarajan (2001). Estimation with Applications to
%       Tracking and Navigation, Section 6.2.
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
h = 0.05; mu = deg2rad(63.4468); T_acc = 20; T_ars = 40;
da = exp(-h/T_acc); dg = exp(-h/T_ars);
ba = [0.1; -0.2; 0.05]; bg = [0.01; 0.02; -0.015];
f = [0.3; -0.2; -9.7]; w = [0.01; -0.02; 0.03];
BIG = 1e12;

% ---- ins_mekf_psi (16 states), no aiding, no compass
Rd = struct('position', BIG*eye(3), 'velocity', BIG*eye(3), 'compass', BIG, ...
    'gravityRefVector', BIG*eye(3), 'pseudoMeasSeaLevel', BIG, ...
    'applicationFlag', false, 'pseudoFlag', false);
x0 = [zeros(6,1); ba; 1; 0; 0; 0; bg]; P0 = eye(15);
Qd = 1e-4*eye(12);
x1 = ins_mekf_psi(x0, P0, mu, h, Qd, Rd, T_acc, T_ars, f, w, []);
nfail = nfail + report('ins_mekf_psi: accelerometer bias', norm(x1(7:9) - da*ba), 1e-12);
nfail = nfail + report('ins_mekf_psi: gyro bias', norm(x1(14:16) - dg*bg), 1e-12);

% ---- ins_mekf (16 states), 6 values: no aiding
Rd2 = Rd; Rd2.magnetometer = BIG*eye(3);
x1 = ins_mekf(x0, P0, mu, h, Qd, Rd2, T_acc, T_ars, [f' w'], [1 0 0]');
nfail = nfail + report('ins_mekf: accelerometer bias', norm(x1(7:9) - da*ba), 1e-12);
nfail = nfail + report('ins_mekf: gyro bias', norm(x1(14:16) - dg*bg), 1e-12);

% ---- ins_euler (15 states), no position aiding; the compass reading has a huge Rd
th0 = [0.1; 0.2; 0.3];
xe = [zeros(6,1); ba; th0; bg]; Pe = eye(15);
Qe = diag([0.01*ones(1,6) 0.1*ones(1,3) 0.001*ones(1,3)]);
Rde = diag([BIG*ones(1,6) BIG]);
x1 = ins_euler(xe, Pe, mu, h, Qe, Rde, T_acc, T_ars, [f; w], th0(3));
nfail = nfail + report('ins_euler: accelerometer bias', norm(x1(7:9) - da*ba), 1e-9);
nfail = nfail + report('ins_euler: gyro bias', norm(x1(13:15) - dg*bg), 1e-9);

% ---- ins_ahrs (9 states), no aiding
xa = [zeros(6,1); ba]; Pa = eye(9);
x1 = ins_ahrs(xa, Pa, mu, h, 1e-3*eye(6), BIG*eye(3), T_acc, f, th0);
nfail = nfail + report('ins_ahrs: accelerometer bias', norm(x1(7:9) - da*ba), 1e-12);

% ---- ins_heave (3 states), no aiding
xh = [0; 0; 0.1];
x1 = ins_heave(xh, eye(3), h, diag([0.01 0.001]), 1, T_acc, mu, 1025, f, 0.1, 0.2);
clear ins_heave
nfail = nfail + report('ins_heave: accelerometer bias', abs(x1(3) - da*0.1), 1e-12);

% ---- quatMEKF, predictor only
q0 = [1; 0; 0; 0];
[~, b1] = quatMEKF(q0, bg, eye(6), h, 1e-3*eye(6), BIG*eye(4), T_ars, [1 0 0]', [f' w']);
nfail = nfail + report('quatMEKF: gyro bias', norm(b1 - dg*bg), 1e-12);

% ---- the decay is the discretisation of db/dt = -b/T: 600 s from the same b
x = x0; for k = 1:round(600/h), x = ins_mekf_psi(x, P0, mu, h, Qd, Rd, T_acc, T_ars, f, w, []); end
nfail = nfail + report('ins_mekf_psi: gyro bias after 600 s against exp(-t/T) b0', norm(x(14:16) - exp(-600/T_ars)*bg), 1e-9);

if nfail > 0
    error('testInsBiasDecay FAILED (%d check(s) wrong)', nfail);
end
disp('testInsBiasDecay PASSED');
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-70s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-70s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
