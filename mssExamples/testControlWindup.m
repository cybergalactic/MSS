function testControlWindup()
% testControlWindup() checks the integral windup of the controllers of the
% GNC folder and of the depth loop of SIMremus100.m when an output limit is
% active.
%
% (1) headingAutopilot, integralSMCheading and PIDnonlinearMIMO have an
%     optional output limit (tauN_max, delta_max, tau_max). With the limit
%     given, the output is limited and the integral state is not updated while
%     the limit holds it (conditional integration). Without it the functions
%     are unchanged (the first check of each function compares with the
%     unlimited call). A caller that only clips the output keeps integrating
%     and overshoots; the closed loops below show the difference.
% (2) SIMremus100.m, control method 1, depth step of 100 m: z_int is frozen
%     while the pitch command is at its limit. The test runs a copy of the
%     script with the plots removed.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   Tanakitkorn (2017), Section 4.1.1, p. 6 (conditional integration).
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
[~,~,Mo] = otter(); M66 = Mo(6,6);

% ---- (1a) headingAutopilot: M r' = tau, step 1.5 rad, fast reference
h = 0.05; wn = 1.5; zeta = 1;
[os0, ~] = headingLoop(M66, wn, zeta, h, []);       % caller clips the moment only
fprintf('  (a caller that only clips the moment overshoots by %.3f rad)\n', os0);
try
    [os1, tauMaxSeen] = headingLoop(M66, wn, zeta, h, 10);
    nfail = nfail + report('headingAutopilot, overshoot with the option (rad)', os1, 0.06);
    nfail = nfail + report('headingAutopilot, output stays within the limit (N m)', max(0, tauMaxSeen - 10), 1e-12);
    clear headingAutopilot
    a = headingAutopilot(0.3, 0, 1, 0.5, M66, wn, zeta, h);
    clear headingAutopilot
    b = headingAutopilot(0.3, 0, 1, 0.5, M66, wn, zeta, h, [], 1e9);
    clear headingAutopilot
    nfail = nfail + report('headingAutopilot, a far limit changes nothing', abs(a - b), 1e-12);
catch ME
    nfail = nfail + failed('headingAutopilot with the limit argument', ME);
end

% ---- (1b) integralSMCheading: Nomoto plant, step 2.5 rad, rudder 20 deg
K = 0.25; T = 1.0; lam = 0.1; phi = 0.1; Kd = 0.5; Ks = 0; dmax = deg2rad(20);
os0 = smcLoop(K, T, lam, phi, Kd, Ks, h, dmax, false);
fprintf('  (a caller that only clips the rudder overshoots by %.3f rad)\n', os0);
try
    os1 = smcLoop(K, T, lam, phi, Kd, Ks, h, dmax, true);
    nfail = nfail + report('integralSMCheading, overshoot with the option (rad)', os1, 0.25);
    clear integralSMCheading
    a = integralSMCheading(0.3, 0, 1, 0, 0, Kd, Ks, lam, phi, K, T, h);
    clear integralSMCheading
    b = integralSMCheading(0.3, 0, 1, 0, 0, Kd, Ks, lam, phi, K, T, h, 1e9);
    clear integralSMCheading
    nfail = nfail + report('integralSMCheading, a far limit changes nothing', abs(a - b), 1e-12);
catch ME
    nfail = nfail + failed('integralSMCheading with the limit argument', ME);
end

% ---- (1c) PIDnonlinearMIMO: Otter 3-DOF, force limited to 30 N, 60 m step
M = Mo([1 2 6],[1 2 6]); D = diag([20 40 30]);
os0 = mimoLoop(M, D, h, []);
fprintf('  (a caller that only clips the force overshoots by %.2f m)\n', os0);
try
    os1 = mimoLoop(M, D, h, [30; 30; 30]);
    nfail = nfail + report('PIDnonlinearMIMO, overshoot with the option (m)', os1, 4.0);
    clear PIDnonlinearMIMO
    eta = [1; 2; 0.1]; nu = [0; 0; 0];
    a = PIDnonlinearMIMO(eta, nu, [3 2 0.1], M, 0.3, 1, 5, h);
    clear PIDnonlinearMIMO
    b = PIDnonlinearMIMO(eta, nu, [3 2 0.1], M, 0.3, 1, 5, h, 1e9);
    clear PIDnonlinearMIMO
    nfail = nfail + report('PIDnonlinearMIMO, a far limit changes nothing', norm(a - b), 1e-12);
catch ME
    nfail = nfail + failed('PIDnonlinearMIMO with the limit argument', ME);
end

% ---- (2) SIMremus100, method 1, depth step 100 m
sim = runRemusCopy(100);
os = max(sim(:,19)) - 100;                          % depth is x(9): column 10 + 9
nfail = nfail + report('SIMremus100, depth overshoot of a 100 m step (m)', os, 2.0);

if nfail > 0
    error('testControlWindup FAILED (%d check(s) wrong)', nfail);
end
disp('testControlWindup PASSED');
end

function [os, tauSeen] = headingLoop(M66, wn, zeta, h, tauMax)
clear headingAutopilot
x = [0; 0]; step = 1.5; tauSeen = 0; N = round(400/h); psiMax = 0;
for k = 1:N
    if isempty(tauMax)
        tau = headingAutopilot(x(1), x(2), step, 0.5, M66, wn, zeta, h, wn);
        tau = max(min(tau, 10), -10);               % the caller clips
    else
        tau = headingAutopilot(x(1), x(2), step, 0.5, M66, wn, zeta, h, wn, tauMax);
    end
    tauSeen = max(tauSeen, abs(tau));
    f = @(s) [s(2); tau/M66];
    k1 = f(x); k2 = f(x + h/2*k1); k3 = f(x + h/2*k2); k4 = f(x + h*k3);
    x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    psiMax = max(psiMax, x(1));
end
os = psiMax - step;
clear headingAutopilot
end

function os = smcLoop(K, T, lam, phi, Kd, Ks, h, dmax, useOption)
clear integralSMCheading
x = [0; 0]; psi_d = 0; r_d = 0; a_d = 0; step = 2.5; psiMax = 0;
for k = 1:round(400/h)
    if useOption
        delta = integralSMCheading(x(1), x(2), psi_d, r_d, a_d, Kd, Ks, lam, phi, K, T, h, dmax);
    else
        delta = integralSMCheading(x(1), x(2), psi_d, r_d, a_d, Kd, Ks, lam, phi, K, T, h);
        delta = max(min(delta, dmax), -dmax);
    end
    f = @(s) [s(2); (K*delta - s(2))/T];
    k1 = f(x); k2 = f(x + h/2*k1); k3 = f(x + h/2*k2); k4 = f(x + h*k3);
    x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    [psi_d, r_d, a_d] = refModel(psi_d, r_d, a_d, step, 1e9, 1.0, 0.5, h, 1);
    psiMax = max(psiMax, x(1));
end
os = psiMax - step;
clear integralSMCheading
end

function os = mimoLoop(M, D, h, tauMax)
clear PIDnonlinearMIMO
eta = [0; 0; 0]; nu = [0; 0; 0]; xMax = 0;
for k = 1:round(600/h)
    if isempty(tauMax)
        tau = PIDnonlinearMIMO(eta, nu, [60 0 0], M, 0.3, 1, 5, h);
        tau = max(min(tau, 30), -30);
    else
        tau = PIDnonlinearMIMO(eta, nu, [60 0 0], M, 0.3, 1, 5, h, tauMax);
    end
    f = @(s) [rotz3(s(3))*s(4:6); M \ (tau - D*s(4:6))];
    s = [eta; nu];
    k1 = f(s); k2 = f(s + h/2*k1); k3 = f(s + h/2*k2); k4 = f(s + h*k3);
    s = s + h/6*(k1 + 2*k2 + 2*k3 + k4);
    eta = s(1:3); nu = s(4:6); xMax = max(xMax, eta(1));
end
os = xMax - 60;
clear PIDnonlinearMIMO
end

function R = rotz3(psi)
R = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
end

function simData = runRemusCopy(zStep)
src = fileread(which('SIMremus100'));
src = rep(src, 'function SIMremus100()', 'function simData = SIMremus100_t()');
src = rep(src, '[ControlFlag, KinematicsFlag] = controlMethod();', 'ControlFlag = 1; KinematicsFlag = 1;');
src = rep(src, 'z_step = 30;', sprintf('z_step = %g;', zStep));
src = rep(src, 'clearvars;', '');
src = rep(src, 'close all;', '');
p = strfind(src, '% PLOTS'); p = p(1);          % cut the script at the plots banner
bar = strfind(src(1:p), '% ====='); cut = bar(end);
src = [src(1:cut-1) 'end' newline];
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp, 's'));
fid = fopen(fullfile(tmp, 'SIMremus100_t.m'), 'w'); fwrite(fid, src); fclose(fid);
addpath(tmp); c2 = onCleanup(@() rmpath(tmp));
set(0, 'DefaultFigureVisible', 'off');
simData = SIMremus100_t();
end

function s = rep(s, a, b)
if ~contains(s, a), error('testControlWindup:copy', 'pattern not found: %s', a); end
s = strrep(s, a, b);
end

function bad = failed(name, ME)
fprintf('%-66s FAILED: %s\n', name, ME.message); bad = 1;
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-66s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-66s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
