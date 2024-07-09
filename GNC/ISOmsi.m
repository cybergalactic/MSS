function [a_z,w_e] = ISOmsi(t)
% ISO 2631-3, 1997 Motion Sickness Index (MSI)
% [a_z, w_e] = ISOmsi(t) computes a_z and w_e as a function of time t.
% 
% INPUTS:
%   t   : Time in hours
%
% OUTPUTS:
%   a_z : Mean of absolute vertical acceleration (m/s^2), i.e. 
%         a_z = mean( abs(a_measured(:,1)) )
%   w_e : Wave encounter frequency (rad/s)
%
% a_z plotted against w_e represents a MSI of 10%, that is 10% sick people
% during t hours.
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
% Date:      5th November 2001
% Revisions: 

% Initialize the index range
f_values = 0.1:0.03:0.63;
num_values = length(f_values);

% Preallocate arrays for efficiency
w_e = zeros(1, num_values);
a_z = zeros(1, num_values);

% Loop over the range of f
for i = 1:num_values
    f = f_values(i);
    
    % Calculate w_e(i)
    w_e(i) = f * 2 * pi;
    
    % Check the condition and assign a_z(i) accordingly
    if (0.1 <= f && f < 0.315)
        a_z(i) = 0.5 * sqrt(2 / t);
    else
        a_z(i) = 0.5 * sqrt(2 / t) * 6.8837 * f^1.67;
    end
end
