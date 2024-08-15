function [LOSangle, LOSrate] = ...
    LOSobserver(LOSangle, LOSrate, LOScommand, h, K_f)
% LOSobserver is compatible with MATLAB and GNU Octave (www.octave.org).
% This function estimates the desired Line-Of-Sight (LOS) angle and
% rate from a discrete-time LOS guidance law command, LOScommand[k], which
% can be computed using the guidance laws LOSchi.m, LOSpsi.m, ILOSpsi.m, 
% ALOSpsi.m, etc. The observer propagates the estimate of the LOS angle 
% according to
%
%   LOSangle[k+1] = LOSangle[k] + h * ( LOSrate[k] + ...
%     K_f * ssa( LOScommand[k] - LOSangle[k]) )
%
% where the LOS yaw rate estimate, LOSrate[k], is computed using numerical
% differentiation, 
% 
%   LOSrate = T_f * s / (T_f * s + 1) * LOSangle
%
% where T_f is the differentiator time constant, which can be determined by
% pole-placement and inspection of the closed-loop system
% 
%    LOSangle/LOScommand = w_n^2 * (T_f*s + 1) / (s^2 + 2*w_n*s + w_n^2)
% 
% If K_f > 0, it follows that T_f = 1 / ( K_f + 2*sqrt(K_f) + 1 ) and that 
% the natural frequency is w_n = K_f + sqrt(K_f). Exact discretization of 
% the observer gives (Fossen 2021, Eqs. B.46-B.47) 
%
%   LOSrate[k] = LOSangle[k][k] - xi[k]
%   xi[k+1] = exp(-h/T_f) * xi[k] + (1 - exp(-h/T_f)) * LOSangle[k]
%
% Inputs:
%   LOSangle:    Estimate of the desired LOS angle at time k
%   LOSrate:     Estimate of the desired LOS angular rate t time k
%   LOScommand:  Commanded LOS angle t time k, computed using ALOSpsi.m,
%                ILOSpsi.m, LOSchi.m, etc.
%   h:           Sampling time (s)
%   K_f:         Observer gain, LOSangle (typically 0.1-0.5)
%
% Outputs:
%   LOSangle:    Updated estimate of the desired LOS angle at time k+1
%   LOSrate:     Updated estimate of the desired LOS rate at time k+1
%
% Calls: ssa.m
%  
% Author:    Thor I. Fossen
% Date:      2024-04-01
%   2024-08-15 : Added pole-placement algorithm for T_f.

% Differentiator time constant
T_f = 1 / ( K_f + 2*sqrt(K_f) + 1 );  

% Internal differentiator state: xi[k]
xi = LOSangle - LOSrate;

% Observer for the LOS angle: LOSangle[k+1]
LOSangle = LOSangle + h * (LOSrate + K_f * ssa( LOScommand - LOSangle) );

% Propagation of the differentiator state: xi[k+1]
PHI = exp(-h/T_f);
xi = PHI * xi + (1 - PHI) * LOSangle;

% LOS rate estimate: LOSrate[k+1]
LOSrate = LOSangle - xi;

end
