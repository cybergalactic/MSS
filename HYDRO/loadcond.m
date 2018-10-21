function loadcond(vessel)
% loadcond    Plots the vessel load condition data
%
%     >> loadcond(vessel)
%
%  Inputs:
%     vessel    : MSS vessel structure
%
% Author:    Thor I. Fossen
% Date:      2005-09-26
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
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

m   = vessel.main.m;
g = 9.81;
Ix = vessel.MRB(4,4);

beta = vessel.main.k44/vessel.main.B;

per = DPperiods(vessel);
T4 = per(4);
w4 = 2*pi/T4;

A44 = interp1(vessel.freqs,reshape(vessel.A(4,4,:),1,length(vessel.freqs)),w4);
M44 = Ix + A44;

GM = 0.1:0.1:15;
T_roll  = 2*pi*sqrt(M44./(m*g*GM));
T_design = 2*pi*sqrt(M44./(m*g*vessel.main.GM_T));

figure(200)
plot(GM,T_roll,'-','linewidth',2)
hold on
plot(vessel.main.GM_T,T4,'r*','MarkerSize',7,'linewidth',2)
plot([0 15],[T_design T_design],'k:')
hold off

legend(sprintf('R_{44} = %3.2f*B',beta),...
    sprintf('Design: GM_T = %3.2f m, T_4 = %3.2f s',vessel.main.GM_T,T_design))
xlabel('GM_T (m)')
ylabel('T_4 (s)')
title('Roll period T_4 (s) as a function of GM_T (m)')