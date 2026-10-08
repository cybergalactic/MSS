function testRemus100Density()
% testRemus100Density() checks that the REMUS 100 model uses one density of
% water. A fully submerged vehicle whose mass is computed from the water it
% displaces (spheroid.m: m = 4/3 * pi * rho * a * b^2, and W = B = m * g)
% has equations of motion that do not depend on the value of rho: the
% rigid-body and added mass, the Coriolis and damping matrices, the
% restoring forces, the propeller, the fins, the lift-drag and the
% cross-flow forces are all proportional to rho, so the accelerations
% nudot = M \ tau do not change when rho changes.
%
% The test evaluates remus100.m and a copy of it with rho = 1025 in place
% of rho = 1026 on 200 seeded states and inputs, and requires the two state
% derivatives to agree to a relative tolerance of 1e-12. It fails if any
% force or mass term inside the model uses a density other than rho.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Author:    Enio Krizman
% Date:      7 October 2026

src = fileread(which('remus100'));
[tok, n] = regexp(src, '^rho = 1026;', 'match', 'lineanchors');
if n == 0 || numel(tok) ~= 1
    error('testRemus100Density: expected one line "rho = 1026;" in remus100.m');
end
src = regexprep(src, '^rho = 1026;', 'rho = 1025;', 'lineanchors');
src = strrep(src, 'function [xdot,U,M] = remus100(', ...
    'function [xdot,U,M] = remus100rho1025(');
folder = tempname; mkdir(folder);
fid = fopen(fullfile(folder, 'remus100rho1025.m'), 'w');
fwrite(fid, src); fclose(fid);
addpath(folder); rehash;
cleanup = onCleanup(@() removeFolder(folder));

rng(1);
maxError = 0;
for k = 1:200
    x = [randn(6,1) .* [1; 0.3; 0.3; 0.3; 0.3; 0.3]
         randn(6,1) .* [10; 10; 5; 0.3; 0.3; pi]];
    ui = [deg2rad(25) * (2*rand-1); deg2rad(25) * (2*rand-1); ...
          1600 * (2*rand-1)];
    Vc = 0.5 * rand; betaVc = 2 * pi * rand; w_c = 0.1 * randn;
    xdot1026 = remus100(x, ui, Vc, betaVc, w_c);
    xdot1025 = remus100rho1025(x, ui, Vc, betaVc, w_c);
    maxError = max(maxError, max(abs(xdot1025(1:6) - xdot1026(1:6))) / ...
        max(abs(xdot1026(1:6))));
end

if maxError < 1e-12
    fprintf('testRemus100Density PASSED: max relative difference %.2e\n', maxError);
else
    error('testRemus100Density FAILED: max relative difference %.2e (tolerance 1e-12)', maxError);
end

end

function removeFolder(folder)
% Remove the temporary copy from the path, then delete it.
rmpath(folder);
rmdir(folder, 's');
end
