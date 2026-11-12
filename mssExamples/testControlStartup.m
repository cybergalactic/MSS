function testControlStartup()
% testControlStartup() checks the start-up of the controllers that keep a
% reference state between calls (persistent variables):
%
%   PIDnonlinearMIMO   the setpoint filter eta_d and
%   headingAutopilot   the reference model state psi_d
%
% started at the origin, not at the measurement. A vehicle that sits on its
% command at (100, 50) m with heading 0.5 rad therefore received a force of
% about -1100 N (the 100 m "error" times Kp) in the first sample, and a yaw
% moment of 82 N m for a heading of 1.14 rad. With the states started at the
% measurement, a vehicle on its command gets zero output.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eq. 12.111-12.112, p. 375 (nonlinear PID with the
%       filtered setpoint).
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
[~,~,M6] = otter();
M = M6([1 2 6],[1 2 6]);
wn = 0.3; zeta = 1; T_f = 5; h = 0.05;

% ---- PIDnonlinearMIMO, 3 DOF, vehicle on its command away from the origin
clear PIDnonlinearMIMO
eta = [100; 50; 0.5]; nu = [0; 0; 0];
tau1 = PIDnonlinearMIMO(eta, nu, eta, M, wn, zeta, T_f, h);
nfail = nfail + report('3 DOF, first force on the command, |tau| (N, N m)', norm(tau1), 1e-9);
mx = 0;
for k = 1:200
    tau = PIDnonlinearMIMO(eta, nu, eta, M, wn, zeta, T_f, h);
    mx = max(mx, norm(tau));
end
nfail = nfail + report('3 DOF, force stays zero for 200 samples, max |tau|', mx, 1e-9);

% ---- 6 DOF
clear PIDnonlinearMIMO
eta6 = [100; 50; 0; 0; 0; 0.5]; nu6 = zeros(6, 1);
tau6 = PIDnonlinearMIMO(eta6, nu6, [100 50 0.5], M6, wn, zeta, T_f, h);
nfail = nfail + report('6 DOF, first force on the command, |tau|', norm(tau6), 1e-9);

% ---- a setpoint step still moves the vehicle towards the new command
clear PIDnonlinearMIMO
t1 = PIDnonlinearMIMO(eta, nu, eta, M, wn, zeta, T_f, h);
t2 = PIDnonlinearMIMO(eta, nu, eta + [10; 0; 0], M, wn, zeta, T_f, h);
for k = 1:20, t2 = PIDnonlinearMIMO(eta, nu, eta + [10; 0; 0], M, wn, zeta, T_f, h); end
nfail = nfail + report('a step of +10 m in x gives a positive surge force (0 = yes)', double(~(t2(1) > 0)), 0.5);
clear PIDnonlinearMIMO

% ---- headingAutopilot
clear headingAutopilot
tauN = headingAutopilot(1.0, 0, 1.0, deg2rad(10), M(3,3), 1.5, 1, h);
nfail = nfail + report('heading autopilot, first moment on the command, |tauN| (N m)', abs(tauN), 1e-9);
mx = 0;
for k = 1:200
    tauN = headingAutopilot(1.0, 0, 1.0, deg2rad(10), M(3,3), 1.5, 1, h);
    mx = max(mx, abs(tauN));
end
nfail = nfail + report('heading autopilot, moment stays zero for 200 samples', mx, 1e-9);
clear headingAutopilot
for k = 1:100, tauN = headingAutopilot(1.0, 0, 2.0, deg2rad(10), M(3,3), 1.5, 1, h); end
nfail = nfail + report('a +1 rad heading step gives a positive moment (0 = yes)', double(~(tauN > 0)), 0.5);
clear headingAutopilot

if nfail > 0
    error('testControlStartup FAILED (%d check(s) wrong)', nfail);
end
disp('testControlStartup PASSED');
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-66s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-66s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
