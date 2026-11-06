function tauN = headingAutopilot(psi, r, psi_ref, r_max, M66, wn, zeta, h, wn_d, tauN_max)
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
% Optional input:
%   wn_d    : Reference-model natural frequency (rad/s), default wn/20
%
%   tauN_max: Yaw moment limit (Nm). When given (use [] for the default wn_d),
%             the output is limited to +-tauN_max and the integral state is
%             not updated while the limit is active (conditional integration,
%             no integral windup). Without it, the output is not limited.
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
%   2026-11-06: Optional yaw moment limit tauN_max with conditional
%               integration; the default is unchanged.

persistent z_int;          % Integral state
persistent a_d r_d psi_d;  % Reference model states

% Initialization of integral state z_int and reference model states a_d, r_d, psi_d
if isempty(z_int)
    z_int = 0;
    a_d = 0;
    r_d = 0;
    psi_d = 0;
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

tauN = -Kp * ssa(psi - psi_d) - Kd * (r - r_d) - Ki * z_int;

% Optional output limit. The integral state is frozen while the limit holds
% the output (conditional integration), otherwise it winds up.
limited = false;
if nargin > 9 && ~isempty(tauN_max)
    if abs(tauN) > tauN_max
        limited = true;
        tauN = sign(tauN) * tauN_max;
    end
end

if ~limited
    z_int = z_int + h * ssa(psi - psi_d);  % Integral state: z_int[k+1]
end

% Propagation of reference model states: psi_d[k+1], r_d[k+1], a_d[k+1]
psi_d = psi_d + h * r_d;
r_d = r_d + h * a_d;
a_d = a_d + h * a_d_dot;

% Maximum desired turning rate
if abs(r_d) > r_max
    r_d = sign(r_d) * r_max;
end

end
