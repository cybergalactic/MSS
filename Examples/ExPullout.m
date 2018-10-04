% ExPullout    Performes a pullout maneuver for two different ships
%
% Author:   Thor I. Fossen
% Date:     25th July 2001
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org

delta_c = 20*pi/180;     % rudder angle for manuver (rad)
h = 0.1;                 % sampling time (sec)

disp('Pullout maneuver for the Mariner class vessel (stable ship)')

% Mariner class cargo ship, cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);        % x  = [ u v r x y psi delta ]' (initial values)
ui = delta_c;           % ui = delta_c 
 
[t,r1,r2] = pullout('mariner',x,ui,h);

% The Esso Osaka tanker (see tanker.m)
disp(' '); disp('Press a key to simulate: The Esso Osaka tanker (unstable ship)'); pause

n = 80;
U = 8.23;
x = [ U 0 0 0 0 0 0 n ]';     % x = [ u v r x y psi delta n ]' (initial values)
n_c     = 80;                 % n_c = propeller revolution in rpm
depth   = 200;                % water depth
ui = [delta_c, n_c, depth];
[t,r1,r2] = pullout('tanker',x,ui,h);
