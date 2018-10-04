% ExTurnCircle  Generates the turning circle for two different ships
%
% Author:   Thor I. Fossen
% Date:     21th July 2001
% Revisions: 11th June 2003 (Thor I. Fossen) rudder execute is changed from 15 deg to 35 deg (IMO regulations)
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

t_final = 700;          % final simulation time (sec)
t_rudderexecute = 100;   % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

disp('Turning cirlce for the Mariner class vessel')

% Mariner class cargo ship, cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = -35*pi/180;   % delta_c = -delta_R at time t = t_rudderexecute
 
[t,u,v,r,x,y,psi,U] = turncircle('mariner',x,ui,t_final,t_rudderexecute,h);

% Container ship (see container.m)
disp(' '); disp('Press a key to simulate: Container Ship'); pause

x     = [8.0 0 0 0 0 0 0 0 0 70]';     % x = [ u v r x y psi delta n ]' (initial values)
delta_c = -35*pi/180;                  % delta_c = -delta_R at time t = t_rudderexecute
n_c     = 80;                          % n_c = propeller revolution in rpm

ui = [delta_c, n_c];
[t,u,v,r,x,y,psi,U] = turncircle('container',x,ui,t_final,t_rudderexecute,h);
