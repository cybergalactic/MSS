function tauN = headingAutopilot(psi, r, psi_ref, r_max, M66, wn, zeta, h)
% tauN = headingAutopilot(psi,r,psi_ref,r_max,M66,wn,zeta,h)
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
% Outputs:
%   tauN    : Autopilot yaw moment command (Nm)
%
% References:
%
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%      Control, 3rd edition, John Wiley & Sons. Ltd., Chichester, UK.

persistent z_int;          % Integral state
persistent a_d r_d psi_d;  % Reference model states

% Initialization of integral state z_int and reference model states a_d, r_d, psi_d
if isempty(z_int)
    z_int = 0;
    a_d = 0;
    r_d = 0;
    psi_d = 0;
end

% Desired jerk (Fossen 2027, Equation 12.10)
wn_d = wn / 20; % Reference signal is 20 times slower than wn
a_d_dot = -wn_d^3 * ssa(psi_d - psi_ref) - 3 * wn_d^2 * r_d ...
    - 3 * wn_d * a_d;

% SISO pole-placement algorithm (Fossen 2027, Algorithm 15.1)
Kp = M66 * wn^2;
Kd = M66 * 2 * zeta * wn;
Ki = 1/10 * Kp * wn;

tauN = -Kp * ssa(psi - psi_d) - Kd * (r - r_d) - Ki * z_int;

z_int = z_int + h * ssa(psi - psi_d);  % Integral state: z_int[k+1]

% Propagation of reference model states: psi_d[k+1], r_d[k+1], a_d[k+1]
psi_d = psi_d + h * r_d;
r_d = r_d + h * a_d;
a_d = a_d + h * a_d_dot;

% Maximum desired turning rate
if abs(r_d) > r_max
    r_d = sign(r_d) * r_max;
end

end
