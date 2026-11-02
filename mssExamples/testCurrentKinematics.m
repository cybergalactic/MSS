function testCurrentKinematics()
% testCurrentKinematics() checks the ocean-current kinematics of the 6-DOF
% models remus100, npsauv, otter, osv and hydroVessel against the equations
% for an irrotational current that is constant in NED (Fossen 2011, Eqs.
% 8.137-8.141 and 8.157):
%
%   v_c^b = R(Theta)' * v_c^n,       nu_c = [v_c^b; 0; 0; 0]
%   d/dt v_c^b = -S(w) * v_c^b,      nu_c_dot = [-S(w) * v_c^b; 0; 0; 0]
%
% where v_c^n = [Vc*cos(betaVc), Vc*sin(betaVc), w_c]' is the current in
% NED, R(Theta) is the rotation matrix from BODY to NED for the full
% attitude and w = [p, q, r]' is the angular velocity.
%
% The test does not look inside the models. In a uniform current that is
% constant in NED, every hydrodynamic force depends on the velocity relative
% to the water, nu_r = nu - nu_c, so the motion relative to the water is
% the motion in still water (Fossen 2011, Eqs. 8.142-8.147). For a state
% x_r = [nu_r; eta] in still water and the same state moving with the
% current, x = [nu_r + nu_c; eta], the two state derivatives must satisfy
%
%   xdot(1:6) - xdot_r(1:6) = nu_c_dot,   attitude rates equal.
%
% The test evaluates this on 200 seeded attitudes, angular velocities,
% velocities, inputs and currents for each model, with the models called
% as they are by default, and requires the error to be below 1e-10
% relative to the still-water accelerations. Two models have one term
% that does not use the relative velocity, and are tested where that term
% has no effect:
%   npsauv: the propeller advance ratio uses the speed over ground u, so
%           the propeller is stopped (n_p = 0).
%   osv:    the added-mass Coriolis matrix is computed from nu, so the
%           vessel drifts with the current (nu_r = 0, w = 0). This checks
%           nu_c but not nu_c_dot.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Section 8.3, pp. 221-225, and Eq. 2.18, p. 22.
%
% Author:    Enio Krizman
% Date:      8 October 2026

N = 200;
tol = 1e-10;
rng(1);

% hydroVessel with the LAUV_marie vessel data
cfg = hydroVesselConfig('LAUV_marie.mat');
d = load(which('LAUV_marie.mat'), 'vessel'); vessel = d.vessel;
vessel.matFile = 'LAUV_marie.mat';
vessel.kappa_4 = cfg.damping.nonlinear_456(1);
vessel.kappa_5 = cfg.damping.nonlinear_456(2);
vessel.kappa_6 = cfg.damping.nonlinear_456(3);
vessel = computeManeuveringModel(vessel, cfg.environment.w0, ...
    cfg.damping.T_1236, cfg.damping.delta_zeta_45, 0);

models = {'remus100', 'remus100 (unit quaternion)', ...
    'npsauv (propeller stopped)', 'otter', ...
    'osv (drifting with the current)', 'hydroVessel'};
failed = false;

for i = 1:numel(models)
    name = models{i};
    maxError = 0; maxAttitude = 0;
    for k = 1:N
        [x_r, call] = sampleCase(name, vessel);
        Theta = x_r(call.attitudeRows);
        w = x_r(4:6);
        R = rotationMatrix(Theta);
        v_cb = R' * call.v_cn;
        nu_c_dot = [-cross(w, v_cb); 0; 0; 0];

        x = x_r; x(1:3) = x(1:3) + v_cb;
        xdot_r = call.still(x_r);
        xdot = call.current(x);

        e = max(abs(xdot(1:6) - xdot_r(1:6) - nu_c_dot)) ...
            / max(1, max(abs(xdot_r(1:6))));
        maxError = max(maxError, e);
        att = call.attitudeRows;
        maxAttitude = max(maxAttitude, max(abs(xdot(att) - xdot_r(att))));
    end
    if maxError < tol && maxAttitude < tol
        status = 'PASSED';
    else
        status = 'FAILED'; failed = true;
    end
    fprintf('%-32s %s: max error %.2e (accelerations), %.2e (attitude rates)\n', ...
        name, status, maxError, maxAttitude);
end

if failed
    error('testCurrentKinematics FAILED (tolerance %.0e)', tol);
else
    fprintf('testCurrentKinematics PASSED (tolerance %.0e)\n', tol);
end

end

function [x, call] = sampleCase(name, vessel)
% One seeded state, input and current, and the two calls of the model.
Vc = rand; betaVc = 2 * pi * (rand - 0.5); w_c = 0.2 * (2*rand - 1);
nu = [1.5 * rand; 0.3 * randn(2,1); 0.3 * randn(3,1)];
pos = [100 * randn(2,1); 20 * rand];
call.attitudeRows = 10:12;

switch name
    case {'remus100', 'remus100 (unit quaternion)'}
        Theta = [deg2rad(30) * (2*rand(2,1) - 1); pi * (2*rand - 1)];
        x = [nu; pos; Theta];
        if strcmp(name, 'remus100 (unit quaternion)')
            x = [nu; pos; euler2q(Theta(1), Theta(2), Theta(3))];
            call.attitudeRows = 10:13;
        end
        ui = [deg2rad(20) * (2*rand(2,1) - 1); 1525 * rand];
        call.still = @(x) remus100(x, ui);
        call.current = @(x) remus100(x, ui, Vc, betaVc, w_c);
    case 'npsauv (propeller stopped)'
        Theta = [deg2rad(30) * (2*rand(2,1) - 1); pi * (2*rand - 1)];
        delta = deg2rad(20) * (2*rand(4,1) - 1);
        x = [nu; pos; Theta; delta; 0];
        ui = [deg2rad(20) * (2*rand(4,1) - 1); 0];
        call.still = @(x) npsauv(x, ui);
        call.current = @(x) npsauv(x, ui, Vc, betaVc, w_c);
    case 'otter'
        w_c = 0;
        Theta = [deg2rad(15) * (2*rand(2,1) - 1); pi * (2*rand - 1)];
        pos(3) = 0.1 * randn;
        x = [nu; pos; Theta];
        n = 100 * (2*rand(2,1) - 1); mp = 25; rp = [0.05 0 -0.35]';
        call.still = @(x) otter(x, n, mp, rp, 0, 0);
        call.current = @(x) otter(x, n, mp, rp, Vc, betaVc);
    case 'osv (drifting with the current)'
        w_c = 0;
        Theta = [deg2rad(15) * (2*rand(2,1) - 1); pi * (2*rand - 1)];
        pos(3) = 0.1 * randn;
        x = [zeros(6,1); pos; Theta];
        ui = [2*rand(4,1) - 1; pi * (2*rand(2,1) - 1)];
        call.still = @(x) osv(x, ui, 0, 0);
        call.current = @(x) osv(x, ui, Vc, betaVc);
    case 'hydroVessel'
        w_c = 0;
        Theta = [deg2rad(30) * (2*rand(2,1) - 1); pi * (2*rand - 1)];
        x = [nu; pos; Theta];
        tau = 20 * randn(6,1);
        call.still = @(x) hydroVessel(x, tau, vessel, 0, 0);
        call.current = @(x) hydroVessel(x, tau, vessel, Vc, betaVc);
end
call.v_cn = [Vc * cos(betaVc); Vc * sin(betaVc); w_c];
end

function R = rotationMatrix(Theta)
% Rotation matrix from BODY to NED: Euler angles, zyx convention (Fossen
% 2011, Eq. 2.18), or unit quaternion q = [eta eps1 eps2 eps3]' (Eq. 2.56).
if length(Theta) == 3
    cphi = cos(Theta(1)); sphi = sin(Theta(1));
    cth  = cos(Theta(2)); sth  = sin(Theta(2));
    cpsi = cos(Theta(3)); spsi = sin(Theta(3));
    R = [ cpsi*cth  -spsi*cphi+cpsi*sth*sphi   spsi*sphi+cpsi*cphi*sth
          spsi*cth   cpsi*cphi+sphi*sth*spsi  -cpsi*sphi+sth*spsi*cphi
          -sth       cth*sphi                  cth*cphi ];
else
    q = Theta / norm(Theta);
    e = q(2:4);
    S = [0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
    R = eye(3) + 2 * q(1) * S + 2 * S^2;
end
end
