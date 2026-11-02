function testGuidanceControlFixes()
% testGuidanceControlFixes() checks crosstrackWpt3D.m and PIDnonlinearMIMO.m.
%
% crosstrackWpt3D: the tracking-error vector p_e = R_y' * R_z' * (p - p1) is
% the vector p - p1 expressed in a rotated frame (the path-tangential frame
% P, as the header of the function states), so it must have the
% same length as p - p1, the waypoint p2 must lie on the x axis of P at the
% distance |p2 - p1|, and a craft displaced from the path to starboard must
% get a positive cross-track error. This is tested on 200 seeded paths and
% positions. A matrix R_z whose (2,1) element has the wrong sign is not a
% rotation (determinant below 1) and fails all three.
%
% PIDnonlinearMIMO: the proportional term uses the yaw error wrapped to
% [-pi, pi) (e(3) = ssa(e(3))). The same must hold for the integral state:
% the heading 3.1 rad with the setpoint -3.1 rad is a 0.083 rad error, not
% a 6.2 rad error. With T_f = h the setpoint filter reaches the setpoint in
% one step, and a setpoint that differs by 2*pi is the same heading, so the
% two control forces must be equal at every step. The yaw force is also
% compared with the control law of the header evaluated with the wrapped
% error.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eq. 2.18, p. 22 (rotation matrices about the
%       principal axes).
%   The control law of PIDnonlinearMIMO.m is the one in its own header.
%
% Author:    Enio Krizman
% Date:      23 October 2026

nfail = 0;

% crosstrackWpt3D
rng(1);
e_norm = 0; e_p2 = 0; neg_cross = 0;
for k = 1:200
    p1 = 100 * randn(3,1);
    p2 = p1 + 100 * randn(3,1);
    p  = p1 + 100 * randn(3,1);
    p_e = crosstrackWpt3D(p2, p1, p);
    e_norm = max(e_norm, abs(norm(p_e) - norm(p - p1)));
    e_p2 = max(e_p2, norm(crosstrackWpt3D(p2, p1, p2) - [norm(p2 - p1); 0; 0]));
end
% Path to the north, a craft 10 m to the east (starboard) of the path, and
% the same geometry turned by other headings: the along-track error is 40 m
% and the cross-track error is +10 m
for psi = [0 30 120 -75] * pi / 180
    R = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    p1 = [5; -3; 2];
    p2 = p1 + R * [100; 0; 0];
    p  = p1 + R * [40; 10; 0];
    p_e = crosstrackWpt3D(p2, p1, p);
    neg_cross = max(neg_cross, norm(p_e - [40; 10; 0]));
end
nfail = nfail + report('crosstrackWpt3D: |p_e| = |p - p1|, max error (m)', e_norm, 1e-10);
nfail = nfail + report('crosstrackWpt3D: p2 on the along-track axis, max error (m)', e_p2, 1e-10);
nfail = nfail + report('crosstrackWpt3D: starboard offset 10 m, max error (m)', neg_cross, 1e-10);

% PIDnonlinearMIMO
M = diag([100 120 50]); wn = 0.5; zeta = 1; h = 0.1; T_f = h;
eta = [0; 0; 3.1]; nu = zeros(3,1);
N = 60;
tau_a = zeros(3,N); tau_b = zeros(3,N); tau_ref = zeros(1,N);
for variant = 1:2
    clear PIDnonlinearMIMO
    eta_ref = [0; 0; -3.1 + (variant - 1) * 2 * pi];
    z3 = 0; eta_d3 = 0;
    for k = 1:N
        tau = PIDnonlinearMIMO(eta, nu, eta_ref, M, wn, zeta, T_f, h);
        if variant == 1
            tau_a(:,k) = tau;
            Kp = M(3,3) * wn^2; Ki = Kp * wn / 10;
            e3 = ssa(eta(3) - eta_d3);
            tau_ref(k) = -(Kp * e3 + Ki * z3);
            z3 = z3 + h * e3;
            eta_d3 = eta_d3 + h * (eta_ref(3) - eta_d3) / T_f;
        else
            tau_b(:,k) = tau;
        end
    end
end
clear PIDnonlinearMIMO
nfail = nfail + report('PIDnonlinearMIMO: setpoints 2*pi apart, max force difference', ...
    max(abs(tau_a(:) - tau_b(:))), 1e-9);
nfail = nfail + report('PIDnonlinearMIMO: yaw force vs the law with wrapped error', ...
    max(abs(tau_a(3,:) - tau_ref)), 1e-9);

if nfail > 0
    error('testGuidanceControlFixes FAILED (%d check(s) wrong)', nfail);
end
disp('testGuidanceControlFixes PASSED');
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-66s PASSED: %.2e\n', name, err); bad = 0;
else
    fprintf('%-66s FAILED: %.2e (tolerance %.0e)\n', name, err, tol); bad = 1;
end
end
