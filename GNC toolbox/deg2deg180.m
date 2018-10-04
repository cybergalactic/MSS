function y = deg2deg180(angle)
% DEG2DEg180 Converts an angle in deg to the interval (-180, 180]
%
% Author:     Roger Skjetne
% Date:       2003-09-05
% Revisions:  2005-04-13 Roger Skjetne  - changed the function s in order
%                                         to account for correct mapping to
%                                         (-180 180].
%
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


% convert angle in deg to the interval <-180, 180)
r = rem(angle+sign(angle)*180,360);
s = sign(sign(angle) + 2*(sign(abs(rem(angle+180,360)/(360)))-1));
% s = sign(angle);

y = r - s*180;