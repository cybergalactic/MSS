function [LOSangle, LOSrate] = LOSobserver(LOSangle, LOSrate, LOScommand, h, K_f)
% Second-order observer for the desired Line-Of-Sight (LOS) angle and rate.
% The command can be computed using LOSchi.m, LOSpsi.m, ILOSpsi.m,
% ALOSpsi.m, or ALOS3D.m. Steps caused by waypoint switching are filtered by
%
%   LOSangle_dot = LOSrate
%   LOSrate_dot  = w_n^2 * ssa(LOScommand - LOSangle) - 2*w_n*LOSrate
%
% where w_n = K_f + sqrt(K_f), retaining the previous natural frequency.
% For a constant command on the selected angle branch, the transfer function
%
%   LOSangle/LOScommand = w_n^2 / (s^2 + 2*w_n*s + w_n^2)
%
% is critically damped. Forward Euler integration advances the angle using
% the input rate. If the caller limits LOSrate between calls, this also
% limits the angle increment. Choose h*w_n < 1 for nonoscillatory response.
% The angle is left unwrapped; ssa selects the shortest command error.
%
% Inputs:
%   LOSangle:    Estimate of the desired LOS angle at time k
%   LOSrate:     Estimate of the desired LOS angular rate at time k
%   LOScommand:  Commanded LOS angle at time k, computed using ALOSpsi.m,
%                ILOSpsi.m, LOSchi.m, etc.
%   h:           Sampling time (s)
%   K_f:         Positive tuning gain (typically 0.1-0.5)
%
% Outputs:
%   LOSangle:    Updated estimate of the desired LOS angle at time k+1
%   LOSrate:     Updated estimate of the desired LOS rate at time k+1
%
% Calls: ssa.m
%
% Author:    Thor I. Fossen
% Date:      2024-04-01
% Revisions:
%   2024-08-15 : Added pole-placement algorithm for T_f.
%   2026-06-27 : Corrected differentiator update to use LOSangle[k].
%   2026-11-03 : Use consistent second-order angle and rate states, removing
%                direct angle-error injection from the angle update.

% -------------------------------------------------------------------------
% Critically damped LOS observer
% -------------------------------------------------------------------------
w_n = K_f + sqrt(K_f);
angleError = ssa(LOScommand - LOSangle);
rate_dot = w_n^2 * angleError - 2*w_n * LOSrate;

% -------------------------------------------------------------------------
% State propagation from k to k+1
% -------------------------------------------------------------------------
LOSangle = LOSangle + h * LOSrate;
LOSrate = LOSrate + h * rate_dot;

end
