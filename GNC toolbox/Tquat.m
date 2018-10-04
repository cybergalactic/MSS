function Tq = Tquat(q)
% Tq = TQUAT(q) computes the quaternion transformation matrix for attitude
%
% Author:   Thor I. Fossen
% Date:     10th September 2010
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2011 Thor I. Fossen and Tristan Perez
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

eta  = q(1); eps1 = q(2); eps2 = q(3); eps3 = q(4); 
 
Tq = 0.5*[...
   -eps1 -eps2 -eps3        
    eta  -eps3  eps2
    eps3  eta  -eps1
   -eps2  eps1  eta   ];
