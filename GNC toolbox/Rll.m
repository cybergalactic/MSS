function R = Rll(l,mu)
% R = Rll(l,mu) computes the Euler angle rotation matrix Rll 
% between ECEF and NED (l = longitude, mu = lattitude)
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
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


cl  = cos(l);
sl  = sin(l);
cmu = cos(mu);
smu = sin(mu);
 
R = [...
   -cl*smu  -sl   -cl*cmu
   -sl*smu   cl   -sl*cmu
    cmu      0    -smu    ];
