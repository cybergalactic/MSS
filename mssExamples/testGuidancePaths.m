function testGuidancePaths()
% testGuidancePaths() checks the path generators of the GNC folder:
%
%   ALOSpsiHermite  y_e must be the perpendicular distance from the vehicle to
%                   the path, not the distance to the nearest path sample. For a
%                   vehicle exactly on the spline the old value was up to half
%                   the sample spacing, with a sign that follows the curvature.
%   order5          with constraint = 1 the limited derivative keeps its own
%                   sign (it was replaced by +3 times the previous derivative,
%                   so the path left the waypoint away from the next one).
%   order3          the knots of the cubic spline are hit to 1e-8 for 80
%                   waypoints (inv(A)*y gave 1e-5).
%
% The reference values are closed forms: y_e = 0 on the path, y_e = d for a
% point at the distance d from a straight path, the derivative of a path piece
% at its start is the constrained chord, the spline passes the waypoints.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2023). An Adaptive Line-of-Sight (ALOS) Guidance Law for Path
%       Following of Aircraft and Marine Craft. IEEE TCST 31(6), 2887-2894.
%   R. Skjetne (2005). The Maneuvering Problem. PhD thesis, NTNU, Eq. 2.11-2.16,
%       pp. 32-33 (Hermite pieces through the waypoints).
%
% Author:    Enio Krizman
% Date:      6 November 2026

nfail = 0;

% ---- (1) ALOSpsiHermite, vehicle on the spline
wpt.pos.x = [0 100 200 300]; wpt.pos.y = [0 100 0 100];
Umax = 2.0; hh = 0.5;
[w_path, x_path, y_path, dx_path, dy_path, pi_h, pp_x, pp_y, Nf] = hermiteSpline(wpt, Umax, hh);
ws = linspace(3, 420, 277); ye = zeros(size(ws));
for i = 1:numel(ws)
    clear ALOSpsiHermite
    xx = ppval(pp_x, ws(i)); yy = ppval(pp_y, ws(i)); idx = floor(ws(i)) + 1;
    [~, ~, ye(i)] = ALOSpsiHermite(w_path, x_path, y_path, dx_path, dy_path, ...
        pi_h, xx, yy, hh, 10, pp_x, pp_y, idx, Nf);
end
nfail = nfail + report('on the spline, max |y_e| (m)', max(abs(ye)), 5e-3);
nfail = nfail + report('on the spline, mean |y_e| (m)', mean(abs(ye)), 1e-3);

% a point d to one side of a straight path, between two samples
wl.pos.x = [0 50 100 150 200]; wl.pos.y = zeros(1, 5);
[w2, x2, y2, dx2, dy2, pih2, ppx2, ppy2, Nf2] = hermiteSpline(wl, Umax, hh);
d = 3.3; errLine = 0;
for xq = [10.37 61.93 133.31]
    for side = [-1 1]
        clear ALOSpsiHermite
        idx = find(x2 <= xq, 1, 'last');
        [~, ~, yq] = ALOSpsiHermite(w2, x2, y2, dx2, dy2, pih2, xq, side*d, hh, 10, ppx2, ppy2, idx, Nf2);
        errLine = max(errLine, abs(abs(yq) - d));
    end
end
nfail = nfail + report('straight path, |y_e| equals the offset 3.3 m, max error (m)', errLine, 1e-9);

% ---- (2) order5 with the constraint, waypoints that turn back
zx = [0 10 11 5 30 31]; zy = zeros(1, 6);
tab = order5(0, 1, zx, zy)';                 % the table has one column per piece
chord = diff(zx);                                   % [10 1 -6 25 1]
dstart = tab(:, 5);                                 % x'(0) of every piece (coefficient a1)
okSign = all(sign(dstart(2:end)) == sign(chord(2:end)'));
nfail = nfail + report('constrained derivative keeps the sign of its chord (0 = yes)', double(~okSign), 0.5);
nfail = nfail + report('piece 3 leaves waypoint 3 towards waypoint 4: x''(0) (-3 expected)', abs(dstart(3) + 3), 1e-9);
dend = zeros(5, 1);                                 % x'(1) of every piece: 5 a5 + 4 a4 + 3 a3 + 2 a2 + a1
for i = 1:5, dend(i) = polyval(polyder(tab(i, 1:6)), 1); end
nfail = nfail + report('first derivative is continuous at the waypoints, max jump', max(abs(dend(1:4) - dstart(2:5))), 1e-9);
tab0 = order5(0, 0, zx, zy)';
nfail = nfail + report('without the constraint nothing changes: x''(0) piece 1', abs(tab0(1, 5) - tab(1, 5)), 1e-12);

% ---- (3) order3 with many waypoints
rng(11); N = 80;
xv = cumsum(10 + 20*rand(1, N + 1)); yv = cumsum(20*randn(1, N + 1));
t3 = order3(xv, yv, 0, 1, 1, 0, 0, 0)';
err = 0;
for i = 1:N
    err = max(err, abs(polyval(t3(i, 1:4), i - 1) - xv(i)));
    err = max(err, abs(polyval(t3(i, 1:4), i) - xv(i + 1)));
    err = max(err, abs(polyval(t3(i, 5:8), i - 1) - yv(i)));
    err = max(err, abs(polyval(t3(i, 5:8), i) - yv(i + 1)));
end
nfail = nfail + report('order3, 80 pieces, max knot error (m)', err, 1e-8);

if nfail > 0
    error('testGuidancePaths FAILED (%d check(s) wrong)', nfail);
end
disp('testGuidancePaths PASSED');
end

function bad = report(name, err, tol)
if err <= tol
    fprintf('%-72s PASSED: %.3g\n', name, err); bad = 0;
else
    fprintf('%-72s FAILED: %.3g (tolerance %.3g)\n', name, err, tol); bad = 1;
end
end
