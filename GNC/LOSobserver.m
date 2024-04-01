function [LOSangle_hat, LOSrate_hat] = ...
    LOSobserver(LOSangle_hat, LOSrate_hat, LOSangle, h, K_f, T_f)
% LOSobserver estimates the desired Line-Of-Sight (LOS) angle and rate 
% from a LOS guidance law command, LOSangle[k]. The observer is
%
%  LOSangle_hat[k+1] = LOSangle_hat[k] + h * ( LOSrate_hat[k] + ...
%    K_f * ssa( LOSangle[k] - LOSangle_hat[k]) )
%
% where the LOS yaw rate is computedusing numerical differentiation, 
% LOSrate = T_f*s / (T_f*s + 1) * LOSangle. Exact discretization gives
% (Fossen 2021, Eqs. B.46-B.47) 
%
%  LOSrate_hat[k] = LOSangle_hat[k][k] - xi[k]
%  xi[k+1] = exp(-h/T_f) * xi[k] + (1 - exp(-h/T_f)) * LOSangle_hat[k]
%
% Inputs:
%   LOSangle_hat: estimate of the LOS angle at time k
%   LOSrate_hat:  estimate of the LOS angular rate t time k
%   LOSangle:     measured LOS angle t time k
%   h:            sampling time (s)
%   K_f:          observer gain, LOSangle (typically 0.1-0.5)
%   T_f:          (OPTIONALLY) differentiator time constant, 
%                     LOSrate = T_f*s / (T_f*s + 1) * LOSangle
%                 If omitted, the LOSrate loop is chosen 5 x faster than 
%                 the LOSangle loop, 1/T_f = 5 * K_f 
%
% Outputs:
%   LOSangle_hat: Updated estimate of the LOS angle at time k+1
%   LOSrate_hat:  Updated estimate of the LOS rate at time k+1
%
% Calls: ssa.m
%  
% Author:    Thor I. Fossen
% Date:      2024-04-01

if nargin == 5       % yaw rate loop is 5 times faster than yaw angle loop 
    T_f = (1/5) * (1/K_f);  % 1/T_f = 5 * K_f  
end

% Internal differentiator state: xi[k]
xi = LOSangle_hat - LOSrate_hat;

% Observer for the LOS angle: LOSangle_hat[k+1]
LOSangle_hat = LOSangle_hat + h * ( LOSrate_hat + ...
    K_f * ssa( LOSangle - LOSangle_hat) );

% Propagation of the differentiator state: xi[k+1]
PHI = exp(-h/T_f);
xi = PHI * xi + (1 - PHI) * LOSangle_hat;

% LOS rate estimate: LOSrate_hat[k+1]
LOSrate_hat = LOSangle_hat - xi;

end
