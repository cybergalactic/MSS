function testRemus100PropellerRollTorque()
% testRemus100PropellerRollTorque() checks the sign of the propeller roll
% moment in remus100.m.
%
% The water exerts a torque K_prop on a propeller that turns with n_p > 0
% (positive about the x axis, clockwise seen from astern), opposing the
% rotation. The hull receives the reaction, -K_prop, so the roll moment on
% the vehicle is negative for n_p > 0 and positive for n_p < 0. This is the
% sign measured on the REMUS vehicle: Prestero (2001, Sec. 4.6.2, Eq. 4.47,
% Table C.9) gives K_prop = -0.543 Nm and a steady roll offset of -5.3 deg at
% 1500 rpm in straight, level flight. The file had tau(4) = +K_prop/10.
%
% The test checks (1) the roll acceleration caused by the propeller at
% phi = 0, p = 0 (the difference between the run at n and the run at 0 rpm,
% all other inputs and states equal) and (2) the sign of the steady roll
% offset reached in straight flight. The magnitude of the scaled-down
% moment (factor 1/10, revision 2023-10-07) is not changed by the fix: the
% steady offset is -0.963 deg at 1500 rpm, against Prestero's -5.3 deg, and
% is only recorded.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. Prestero (2001). Verification of a six-degree of freedom simulation
%       model for the REMUS autonomous underwater vehicle. MSc thesis,
%       MIT/WHOI. Sec. 4.6.2, Eq. 4.47.
%
% Author:    Enio Krizman
% Date:      2 November 2026

fail = 0;

% (1) Roll acceleration due to the propeller, same state, n versus 0 rpm
x = zeros(12,1); x(1) = 1.5;
xd0 = remus100(x, [0 0 0]');
for n = [1000 1500 1525 -500 -1000]   % n > 0 above the windmilling range of the linear K_Q fit
    xd = remus100(x, [0 0 n]');
    dp = xd(4) - xd0(4);
    ok = sign(dp) == -sign(n);
    fprintf('n = %5d rpm: propeller roll acceleration %+9.5f rad/s^2  (hull torque opposes rotation: %s)\n', ...
        n, dp, passFail(ok));
    fail = fail + ~ok;
end

% (2) Steady roll offset in straight flight
h = 0.02;
for n = [1500 -1000]
    x = zeros(12,1); x(1) = sign(n) * 1.5;
    for k = 1:round(1500/h)
        x = x + h * remus100(x, [0 0 n]');
    end
    phi = rad2deg(x(10));
    ok = sign(phi) == -sign(n);
    fprintf('n = %5d rpm: steady roll offset %+.3f deg  (negative for n > 0: %s)\n', n, phi, passFail(ok));
    fail = fail + ~ok;
    if n == 1500
        fprintf('             recorded: Prestero measured -5.3 deg at 1500 rpm (magnitude not part of the test)\n');
    end
end

if fail > 0
    error('testRemus100PropellerRollTorque: %d check(s) FAILED', fail);
end
disp('testRemus100PropellerRollTorque PASSED');
end

function s = passFail(ok)
if ok, s = 'ok'; else, s = 'FAIL'; end
end
