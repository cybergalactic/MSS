function [a_z,w_e] = ISOmsi(t)
% [a_z,w_e] = ISOmsi(t) computes a_z and w_e as a function of time t
% 
% t  : time in hours
% a_z: mean of absolute vertical acceleration (m/s^2), i.e. a_z = mean(abs(a_measured(:,1)))
% w_e: enconter frequency (rad/s), see encounter.m.
%
% a_z plotted against w_e represents a MSI of 10%, that is 10% sick people
% during t hours.
%
% Refs.  - A. R. J. M. Lloyd (1989). Seakeeping Behaviour in Rough Water. Ellis Horwoowd Ltd.
%        - E. V. Lewis (Ed.) (1989). Principles of Naval Architecture. Vol III Motions in
%             Waves and Controllability, 2nd. ed., SNAME. 
%        - J.F. O'Hanlon and M. E. McCauley (1974). Motion Sickness Incidence as a Function of
%             Vertical Sinusoidal Motion. Aerospace Medicine AM-45(4):366-369.
%
% Author:    Thor I. Fossen
% Date:      5th November 2001
% Revisions: 
%
% ISO 2631-3, 1997 Motion Sickness Index (MSI)
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

i = 1;
for f = 0.1:0.03:0.63,
    w_e(i) = f*2*pi;
    if (0.1 <= f & f < 0.315),
        a_z(i) = 0.5*sqrt(2/t);
    else
        a_z(i) = 0.5*sqrt(2/t)*6.8837*f^1.67;
    end
    i = i+1;
end