function testOtterLegacyTwoPontoonCrossFlow()
% testOtterLegacyTwoPontoonCrossFlow checks the strip-theory cross-flow drag
% of otter_legacy.m against the two-pontoon sum used in otter.m.
%
% otter.m sums the drag of each pontoon, evaluated at its lateral offset
% y = -y_pont and y = +y_pont, and maps each contribution to the center of
% origin with H' (lines 297-304). otter_legacy.m must do the same. This test
% reads the strip block of otter_legacy.m as written in the file, evaluates it
% for several relative velocities, and compares the result with the two-pontoon
% sum built from crossFlowDrag and Hmtrx. A pure sway velocity must give twice
% the single-pontoon sway drag.
%
% The test fails while the block evaluates a single centerline pontoon and
% passes once both pontoons are summed.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'LIBRARY')), fullfile(root,'HYDRO'));

src = fileread(fullfile(root,'CRAFT','USV','models','otter_legacy.m'));
iA = strfind(src,'% Strip theory');
iB = strfind(src,'% Payload expressed in NED');
assert(numel(iA) == 1 && numel(iB) == 1, ...
    'Strip block markers not found exactly once in otter_legacy.m');
block = src(iA:iB-1);

% Hull data from the file itself; the draft is a representative value
L      = readScalar(src,'L');
B_pont = readScalar(src,'B_pont');
y_pont = readScalar(src,'y_pont');
T      = 0.13;

% Relative velocities: surge, sway, heave, roll, pitch, yaw
cases = [ 0.8  0.3  0.1  0.2  0.1  0.4
          0    0.5  0    0    0    0
          0    0    0    0    0    0.3
          0.4 -0.2  0.05 0    0.1 -0.2 ];

nFail = 0;
for k = 1:size(cases,1)
    nu_r = cases(k,:)';

    tau_file = evalStrip(block, L, B_pont, T, y_pont, nu_r);

    tau_ref = zeros(6,1);
    for y = [-y_pont, y_pont]
        H_pont = Hmtrx([0 y 0]);
        tau_ref = tau_ref + H_pont' * crossFlowDrag(L,B_pont,T,H_pont*nu_r);
    end

    err   = norm(tau_file - tau_ref);
    scale = max(norm(tau_ref), eps);
    fprintf('case %d: |tau_file - tau_two_pontoon| / |tau| = %.3e\n', k, err/scale);
    if err > 1e-9 * scale
        nFail = nFail + 1;
    end
end

% Pure sway: each pontoon sees the same sway velocity, so the sway drag is doubled
nu_sway = cases(2,:)';
tau_sway_file = evalStrip(block, L, B_pont, T, y_pont, nu_sway);
tau_single    = crossFlowDrag(L,B_pont,T,nu_sway);
errSway = abs(tau_sway_file(2) - 2*tau_single(2));
fprintf('pure sway: tau_file(2) = %.6g N, 2 x single pontoon = %.6g N\n', ...
    tau_sway_file(2), 2*tau_single(2));
if errSway > 1e-9 * max(abs(tau_single(2)), eps)
    nFail = nFail + 1;
end

assert(nFail == 0, ...
    'otter_legacy.m cross-flow block differs from the two-pontoon sum in %d check(s)', nFail);
fprintf('testOtterLegacyTwoPontoonCrossFlow: all checks passed\n');
end

function tau_crossflow = evalStrip(block, L, B_pont, T, y_pont, nu_r)
% Evaluates the strip block of otter_legacy.m in this workspace
eval(block);
end

function v = readScalar(src, name)
% Reads a scalar assignment "name = value;" from the model source
tok = regexp(src, ['^' name '\s*=\s*([0-9.]+)\s*;'], 'tokens', 'once', 'lineanchors');
assert(~isempty(tok), 'Could not read %s from otter_legacy.m', name);
v = str2double(tok{1});
end
