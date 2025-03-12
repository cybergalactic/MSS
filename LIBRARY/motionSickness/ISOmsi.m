function [a_z,w_e] = ISOmsi(t)
% ISO 2631-3, 1997 Motion Sickness Index (MSI)
% [a_z, w_e] = ISOmsi(t) computes a_z and w_e as a function of time t.
% 
% INPUTS:
%   t   : Time in hours
%
% OUTPUTS:
%   a_z : Mean of absolute vertical acceleration (m/s²)
%   w_e : Wave encounter frequency (rad/s)
%
% a_z plotted against w_e represents an MSI of 10%, meaning 10% of people
% will experience motion sickness after `t` hours of exposure.
%
% Refs.  
%   - A. R. J. M. Lloyd (1989). Seakeeping Behaviour in Rough Water. Ellis
%     Horwoowd Ltd.
%   - E. V. Lewis (Ed.) (1989). Principles of Naval Architecture. Vol III 
%     Motions in Waves and Controllability, 2nd. ed., SNAME. 
%   - J.F. O'Hanlon and M. E. McCauley (1974). Motion Sickness Incidence as 
%     a Function of Vertical Sinusoidal Motion. Aerospace Medicine 
%     AM-45(4):366-369.
%
% Author:    Thor I. Fossen
% Date:      2001-11-5
% Revisions: Optimized for vectorized computation

% Define encounter frequency range
f_values = (0.1:0.03:0.63)'; % Column vector

% Convert to rad/s
w_e = 2 * pi * f_values; % Wave encounter frequency (rad/s)

% Avoid division by zero (in case t is very small)
t = max(t, eps); % Ensures no divide-by-zero errors

% Compute a_z based on the condition
a_z = 0.5 * sqrt(2 / t) * ones(size(f_values)); % Default value
high_freq_idx = f_values >= 0.315; % Indices where f >= 0.315

% Apply the frequency-dependent scaling for f >= 0.315
a_z(high_freq_idx) = a_z(high_freq_idx) .* (6.8837 .* f_values(high_freq_idx).^1.67);

end
