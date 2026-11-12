function tauN = headingAutopilot( ...
    psi, r, psi_ref, r_max, M66, wn, zeta, h, wn_d, tauN_max)
% tauN = headingAutopilot(psi,r,psi_ref,r_max,M66,wn,zeta,h)
% tauN = headingAutopilot(psi,r,psi_ref,r_max,M66,wn,zeta,h,wn_d)
% tauN = headingAutopilot(psi,r,psi_ref,r_max,M66,wn,zeta,h,wn_d,tauN_max)
% SISO PID pole-placement algorithm (Fossen 2027, Algorithm 15.1) for
% heading control.
%
%    d/dt z_int = ssa(psi - psi_d)
%
%    tauN = -( Kp * ssa(psi - psi_d) + Ki * z_int ) - Kd * (r - r_d)
%
%    Kp = M(6,6) * wn * wn,
%    Kd = M(6,6) * 2 * zeta * wn
%    Ki = 1/10 * Kp * wn
%
% Persistent variables:
%    The integral state z_int is a persistent variable that should be cleared by
%    adding:
%
%    clear headingAutopilot;
%
% in the top of the script calling the function.
%
% Inputs:
%   psi     : Yaw angle (rad)
%   r       : Yaw rate (rad/s)
%   psi_ref : Desired yaw angle (rad)
%   r_max   : Maximum desired yaw rate (rad/s)
%   M66     : Moment of inertia including hydrodynamic moment of inertia in yaw
%   wn      : Closed-loop natural frequency in yaw (rad/s)
%   zeta    : Closed-loop relative damping ratio in yaw (-)
%   h       : Sampling time (s)
%
% Optional inputs:
%   wn_d     : Reference-model natural frequency (rad/s), default wn/20.
%   tauN_max : Positive yaw moment limit (Nm), default Inf. A finite limit
%              enables anti-windup. Use [] to select a default.
%
% Outputs:
%   tauN    : Autopilot yaw moment command (Nm)
%
% References:
%
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%      Control, 3rd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-10-25
% Revisions:
%   2026-11-06: The reference model state psi_d starts at the measured heading
%               instead of zero.
%   2026-11-12: Add optional yaw moment saturation and anti-windup.

persistent z_int;          % Integral state
persistent a_d r_d psi_d;  % Reference model states

% Initialization of integral state z_int and reference model states a_d, r_d,
% psi_d. The reference model starts at the measured heading.
if isempty(z_int)
    z_int = 0;
    a_d = 0;
    r_d = 0;
    psi_d = psi;
end

% Select reference bandwidth independently of the PID bandwidth.
if nargin < 9 || isempty(wn_d)
    wn_d = wn / 20;
end

% Desired jerk (Fossen 2027, Equation 12.10)
a_d_dot = -wn_d^3 * ssa(psi_d - psi_ref) - 3 * wn_d^2 * r_d ...
    - 3 * wn_d * a_d;

% SISO pole-placement algorithm (Fossen 2027, Algorithm 15.1)
Kp = M66 * wn^2;
Kd = M66 * 2 * zeta * wn;
Ki = 1/10 * Kp * wn;

e_psi = ssa(psi - psi_d);
tauN = -Kp * e_psi - Kd * (r - r_d) - Ki * z_int;
integral_step = h * e_psi;

if nargin >= 10 && ~isempty(tauN_max) && abs(tauN) >= tauN_max
    % Anti-windup blocks increments of -Ki*z_int that increase saturation.
    if tauN * integral_step < 0
        integral_step = 0;
    end
    tauN = sign(tauN) * tauN_max;
end

z_int = z_int + integral_step;  % Integral state: z_int[k+1]

% Propagation of reference model states: psi_d[k+1], r_d[k+1], a_d[k+1]
psi_d = psi_d + h * r_d;
r_d = r_d + h * a_d;
a_d = a_d + h * a_d_dot;

% Maximum desired turning rate
if abs(r_d) > r_max
    r_d = sign(r_d) * r_max;
end

end
