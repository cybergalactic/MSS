function We = encounter(Ww,U,beta)
% We = ENCOUNTER(Ww,U,beta) computes the encounter frequency We (rad/s) as a function
%       of wave peak frequeny Ww (rad/s), vessel speed U (m/s) and wave direction 
%       beta, 0 (deg) for following seas and 180 (deg) for head seas.
%       
% We = abs(Ww - Ww.^2*U*cos(beta*pi/180)/g);
% size(We)=size(Ww)
%
% Author:   Thor I. Fossen
% Date:     2nd November 2001
% Revisions: 
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


g = 9.8; % acceleration of gravity
We = abs(Ww - Ww.^2*U*cos(beta*pi/180)/g);