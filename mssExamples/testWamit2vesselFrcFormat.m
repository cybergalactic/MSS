function testWamit2vesselFrcFormat(titleLine, mass)
% testWamit2vesselFrcFormat(titleLine, mass) checks that wamit2vessel.m tells the two
% WAMIT *.frc layouts apart by structure, not by the magnitude of any one
% value.
%
% The Force Control File (WAMIT user manual, section 4, "the *.frc input
% file") has two layouts. Alternative 1: IOPTN flags, then VCG, then the
% 3x3 radius-of-gyration matrix (three rows of three numbers). Alternative
% 2: IOPTN flags, then RHO, then the centre of gravity [xg yg zg], an
% IMASS flag, then the 6x6 rigid-body mass matrix (one row of one number,
% the flag, followed by six rows of six numbers). The MSS configuration
% guide (Fossen, "Description of MSS Vessel Models", 19 June 2008, sec.
% 2.1.1) prints both: fpso.frc and the generic myvessel.frc example use
% alternative 1; tanker.frc and semisub.frc use alternative 2.
%
% wamit2vessel.m (within this file: "if frc(5,1) < 1000") told the two
% apart by testing whether the fifth data row's first column -- the
% rigid-body mass in alternative 2 -- is below 1000. MSS's own three
% WAMIT vessels pass this test only because their masses, in the units
% the files use, happen to be in the millions. A lighter body in
% alternative 2 (any vessel whose WAMIT run used a convenient model scale)
% has a mass below 1000 in those same units and is misread as alternative
% 1: the density is forced to 1025 regardless of the file's own RHO, and
% the radius-of-gyration fields are populated from the wrong rows of the
% 6x6 mass matrix. Three vessels in a public modelling suite's fleet hit
% exactly this case: masses of 127.92, 13.94 and 13.15 (their own units),
% each run at a density of 1000 given explicitly in the file.
%
% This test builds a synthetic alternative-2 WAMIT output set for a unit
% body with RHO = 1000 and a mass of 500 (both under the 1000 threshold),
% a 6x6 mass matrix with known, distinct diagonal entries, and the
% minimum *.pot, *.out, *.1, *.3, *.4 and *.9 files wamit2vessel.m needs
% to run to completion (traced by reading the function; no WAMIT
% installation is used). It calls wamit2vessel and checks the resulting
% vessel.main.rho, the radii of gyration, vessel.MRB, and the density-
% scaled vessel.A and vessel.C against the values that follow directly
% from the file's own numbers and the WAMIT manual's scaling formulas
% (added mass/damping: rho * ULEN^k, p. 4-3; restoring: rho * g * ULEN^k,
% p. 4-2, both cited already in wamit2vessel.m).
%
% On the unmodified file this test fails: rho reads back as 1025 and the
% radii of gyration, A and C are built from the wrong density and the
% wrong rows. After detecting the layout from the file's structure (the
% row after VCG/CG has three numbers for alternative 1, one for
% alternative 2), it passes.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% References:
%   WAMIT, Inc. WAMIT User Manual, version 7.x, section 4 (Force Control
%       File, *.frc).
%   T. I. Fossen (2008). Description of MSS Vessel Models: Configuration
%       Guidelines for Hydrodynamic Codes. Section 2.1.1, pp. 5-6
%       (FRC-file alternatives 1 and 2), section 2.2 example files
%       (fpso.frc p. 13, tanker.frc p. 16, semisub.frc pp. 9-10).
%
% Author:    Enio Krizman
% Date:      8 October 2026 (revised 2 November 2026 for MSS d78d97e)

% Defaults: the original case (title text, body mass 500). The blank-title
% case (testWamit2vesselFrcFormatBlankTitle) passes '' and 30.
if nargin < 1, titleLine = 'synthetic unit body FRC file (alternative 2, mass < 1000)'; end
if nargin < 2, mass = 500; end
label = 'testWamit2vesselFrcFormat';
if isempty(strtrim(titleLine)), label = [label ' (blank title line)']; end

tol = 1e-9;
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
name = 'unitbody';
f = fullfile(folder, name);

% ---- synthetic inputs (all self-defined, none from a real vessel) ----
rho = 1000; g = 9.81; ulen = 10;
xg = 0.1; yg = 0; zg = 0.2;      % centre of gravity, out-file convention
I44 = 50; I55 = 80; I66 = 90;    % distinct diagonal rigid-body inertias
MRB_in = diag([mass mass mass I44 I55 I66]);
Abar11 = 2.5; Bbar11 = 0.3;      % added mass / damping, DOF (1,1), all periods
C33 = 100; C44 = 200; C55 = 300; % restoring stiffness diagonal terms
T_draught = 1; Lpp = ulen; Boa = 5;

writeFrc([f '.frc'], titleLine, rho, xg, yg, zg, MRB_in);
writePot([f '.pot']);
writeOut([f '.out'], g, ulen, xg, yg, zg, C33, C44, C55);
writeAddedMass([f '.1'], Abar11, Bbar11);
writeRaoOrWD([f '.3'], false);
writeRaoOrWD([f '.4'], false);
writeRaoOrWD([f '.9'], true);

cwd = pwd;
cleanupCwd = onCleanup(@() cd(cwd)); %#ok<NASGU>
cd(folder);
% wamit2vessel saves <name>.mat beside itself; remove it afterwards
cleanupMat = onCleanup(@() delete(fullfile(fileparts(which('wamit2vessel')), [name '.mat'])));
% damping inputs T_126, delta_zeta_345, nonlinear_456 given, so no prompt
vessel = wamit2vessel(name, T_draught, Lpp, Boa, '0000', [5 5 5], [0 0.1 0], [5 0 0]);

% ---- expected values, derived from the inputs above and the manual's
% scaling formulas, not copied from any run of this test ----
expected.rho = rho;
expected.k44 = sqrt(I44/mass);
expected.k55 = sqrt(I55/mass);
expected.k66 = sqrt(I66/mass);
expected.MRB = MRB_in([1 2 3 4 5 6],[1 2 3 4 5 6]);  % FSD axes, Tscale*M*Tscale
% Tscale = diag([1 -1 -1 1 -1 -1]): a pure diagonal sign congruence leaves
% every diagonal entry, and every (i,i)-(j,j) product-of-equal-sign pair
% on this purely diagonal MRB_in, unchanged.
expected.A11 = Abar11 * rho * ulen^3;      % manual p. 4-3, k=3 for i,j<=3
rho_g = rho * g;
expected.C33 = C33 * rho_g * ulen^2;       % manual p. 4-2
expected.C44 = C44 * rho_g * ulen^4;
expected.C55 = C55 * rho_g * ulen^4;

fails = {};
if abs(vessel.main.rho - expected.rho) > tol
    fails{end+1} = sprintf('rho: got %.6g, expected %.6g', vessel.main.rho, expected.rho);
end
if abs(vessel.main.k44 - expected.k44) > tol
    fails{end+1} = sprintf('k44: got %.6g, expected %.6g', vessel.main.k44, expected.k44);
end
if abs(vessel.main.k55 - expected.k55) > tol
    fails{end+1} = sprintf('k55: got %.6g, expected %.6g', vessel.main.k55, expected.k55);
end
if abs(vessel.main.k66 - expected.k66) > tol
    fails{end+1} = sprintf('k66: got %.6g, expected %.6g', vessel.main.k66, expected.k66);
end
if max(abs(diag(vessel.MRB) - diag(expected.MRB))) > tol
    fails{end+1} = sprintf('MRB diagonal: got %s, expected %s', ...
        mat2str(diag(vessel.MRB)'), mat2str(diag(expected.MRB)'));
end
if any(abs(vessel.A(1,1,:) - expected.A11) > tol * max(1,abs(expected.A11)))
    fails{end+1} = sprintf('A(1,1,:): got %s, expected %.6g', ...
        mat2str(squeeze(vessel.A(1,1,:))'), expected.A11);
end
if any(abs(vessel.C(3,3,:) - expected.C33) > tol * max(1,abs(expected.C33)))
    fails{end+1} = sprintf('C(3,3,:): got %s, expected %.6g', ...
        mat2str(squeeze(vessel.C(3,3,:))'), expected.C33);
end
if any(abs(vessel.C(4,4,:) - expected.C44) > tol * max(1,abs(expected.C44)))
    fails{end+1} = sprintf('C(4,4,:): got %s, expected %.6g', ...
        mat2str(squeeze(vessel.C(4,4,:))'), expected.C44);
end
if any(abs(vessel.C(5,5,:) - expected.C55) > tol * max(1,abs(expected.C55)))
    fails{end+1} = sprintf('C(5,5,:): got %s, expected %.6g', ...
        mat2str(squeeze(vessel.C(5,5,:))'), expected.C55);
end

if isempty(fails)
    fprintf('%s PASSED\n', label);
else
    fprintf('%s FAILED:\n', label);
    for i = 1:numel(fails)
        fprintf('  %s\n', fails{i});
    end
    error('%s FAILED (%d check(s) wrong)', label, numel(fails));
end

end

% =========================================================================
function writeFrc(file, titleLine, rho, xg, yg, zg, MRB)
% Alternative 2: title line (may be blank), IOPTN flags, RHO, [xg yg zg],
% IMASS flag, 6x6 mass matrix, IDAMP flag, ISTIFF flag (no external
% damping/stiffness supplied). Structure per the MSS guide's tanker.frc /
% semisub.frc examples (pp. 9-10, 16) and the WAMIT manual's
% alternative-2 layout.
fid = fopen(file,'w');
fprintf(fid, '%s\n', titleLine);
fprintf(fid, '1 0 1 1 0 0 0 1 0\n');
fprintf(fid, '%.1f\n', rho);
fprintf(fid, '%.4f %.4f %.4f\n', xg, yg, zg);
fprintf(fid, '1\n');
for i = 1:6
    fprintf(fid, '%s\n', sprintf('%.4f ', MRB(i,:)));
end
fprintf(fid, '0\n0\n');
fclose(fid);
end

function writePot(file)
% Only requirement read by wamit2vessel.m: a line that is exactly "19"
% (19 headings 0:10:180, MSS's requirement, guide p. 5).
fid = fopen(file,'w');
fprintf(fid, 'unitbody\n0\n0\n-1 0.00 0.00 0.00 0.00\n0 0 0\n1 1\n');
fprintf(fid, '1 1 1 1 1 1\n2\n-1 0\n19\n');
fprintf(fid, '0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180\n');
fclose(fid);
end

function writeOut(file, g, ulen, xg, yg, zg, C33, C44, C55)
% Mimics the fixed-column *.out fields wamit2vessel.m slices by character
% position (gravity, length scale, volumes, C(3,3..5), centre of buoyancy,
% centre of gravity), copying the column layout of MSS's own fpso.out.
% Four dummy "period" rows separate the POTEN header from the Gravity
% line so Nfreqs (= count - countstart - 2) equals the 3 periods in the
% *.1 file, as it does in a real WAMIT run (one row per period computed).
fid = fopen(file,'w');
fprintf(fid, 'Synthetic unit body for testWamit2vesselFrcFormat\n');
fprintf(fid, 'POTEN run date and starting time:        08-Oct-2026  --  00:00:00\n');
fprintf(fid, '   Period       Time           RAD      DIFF  (max iterations)\n');
fprintf(fid, '   -1.0000    00:00:01          -1\n');
fprintf(fid, '    0.0000    00:00:01          -1\n');
fprintf(fid, '    1.0000    00:00:02          -1      -1\n');
fprintf(fid, '    2.0000    00:00:03          -1      -1\n');
fprintf(fid, ' Gravity:     %.5f                Length scale:          %.5f\n', g, ulen);
fprintf(fid, ' Volumes (VOLX,VOLY,VOLZ):         0.5000       0.5000       0.5000\n');
fprintf(fid, ' Center of Buoyancy (Xb,Yb,Zb):      0.000000     0.000000    -0.300000\n');
fprintf(fid, ' C(3,3),C(3,4),C(3,5):  %.4f       0.0000       0.0000\n', C33);
fprintf(fid, ' C(3,3),C(3,4),C(3,5):  %.4f       0.0000       0.0000\n', C44);
fprintf(fid, ' C(3,3),C(3,4),C(3,5):  %.4f       0.0000\n', C55);
fprintf(fid, ' Center of Gravity  (Xg,Yg,Zg):      %.6f     %.6f     %.6f\n', xg, yg, zg);
fclose(fid);
end

function writeAddedMass(file, Abar11, Bbar11)
% *.1: period i j Abar Bbar (manual p. 4-8). Periods -1 (infinite
% frequency), 0 (zero frequency) and one positive period, all DOF (1,1)
% to keep the test to the minimum the detection fault touches.
fid = fopen(file,'w');
for per = [-1 0 1]
    fprintf(fid, '%d 1 1 %.4f %.4f\n', per, Abar11, Bbar11);
end
fclose(fid);
end

function writeRaoOrWD(file, isWD)
% *.3/*.4 (force/motion RAO, manual p. 4-8, 7 columns) or *.9 (wave drift
% from pressure, 8 columns): all 19 MSS headings, one period, zero
% amplitude -- wamit2vessel.m's 0-180 to 180-360 mirror step indexes
% fixed columns up to 36 regardless of content, so all 19 headings must
% be present for the function to run to completion.
fid = fopen(file,'w');
headings = 0:10:180;
for b = headings
    for dof = 1:6
        if isWD
            fprintf(fid, '1.0 %.1f %.1f %d 0.0 0.0 0.0 0.0\n', b, b, dof);
        else
            fprintf(fid, '1.0 %.1f %d 0.0 0.0 0.0 0.0\n', b, dof);
        end
    end
end
fclose(fid);
end
