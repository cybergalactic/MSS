function testRemusDepthLoop()
% testRemusDepthLoop() checks the outer depth loop of SIMremus100.m (control
% method 1, gradient-based depth control). The loop forms the residual
%
%   sigma = z_dot - alpha_z,   alpha_z = z_d_dot - Kp_z*(e_z + z_int/T_z),
%
% so that sigma = 0 gives e_z' = -Kp_z*(e_z + z_int/T_z) for e_z = z - z_d.
% The script wrote sigma without the derivative z_d_dot of the low-pass
% filtered depth command, which therefore acts as a disturbance of the error
% loop during a depth step (peak error 6.56 m for a 30 m step against 3.58 m
% with the term). The test runs a copy of the script with the plots removed,
% a 30 m step at 200 s, and checks the peak depth error and the overshoot.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% The derivative of the filtered command, z_d' = wn_d_z*(z_ref - z_d), is the
% exact derivative of the discrete first-order filter makeLowPass(wn_d_z, h).
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;
rng(0);
sim = runRemusCopy(30);
h = 0.05; i0 = round(200/h) + 1;
z_d = sim(:,1); zn = sim(:,19);
peakErr = max(abs(zn(i0:end) - z_d(i0:end)));
over = max(zn) - 30;
fprintf('  peak |z - z_d| %.3f m, final z %.3f m\n', peakErr, zn(end));
nfail = nfail + report('peak depth error of a 30 m step (m)', peakErr, 4.5);
nfail = nfail + report('overshoot of the 30 m step (m)', max(over, 0), 0.3);

if nfail > 0
    error('testRemusDepthLoop FAILED (%d check(s) wrong)', nfail);
end
disp('testRemusDepthLoop PASSED');
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

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-66s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-66s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
