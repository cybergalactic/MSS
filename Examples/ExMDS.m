% ExMD       Plots the step response of a 2nd-order mass-damper system
% Author:    Thor I. Fossen
% Date:      16th June 2001
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

wn = 1;

subplot(211)
t = 0:0.01:20;
z = 0.5; sys = tf([wn*wn],[1 2*z*wn wn*wn]); step(sys,t)
hold on
z = 1.0; sys = tf([wn*wn],[1 2*z*wn wn*wn]); step(sys,t)
z = 2.0; sys = tf([wn*wn],[1 2*z*wn wn*wn]); step(sys,t)
hold off
legend('\zeta = 0.5','\zeta = 1.0','\zeta = 2.0')

subplot(212)
t = 0:0.01:50;
z = 0.1; sys = tf([wn*wn],[1 2*z*wn wn*wn]); step(sys,t)
hold on
sys = tf([wn*wn],[1 0 wn*wn]); step(sys,t)
hold off
legend('\zeta = 0.1','\zeta = 0.0')
