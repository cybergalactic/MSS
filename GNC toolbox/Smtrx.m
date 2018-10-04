function S = Smtrx(a)
% S = Smtrx(a) computes the 3x3 vector skew-symmetric matrix S(a) = -S(a)'.
% The corss product satisfies: a x b = S(a)b. The inverse can be computed
% as: vex(S(a)) = a 
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 6th August 2011 - updated documentation
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
 
S = [    0  -a(3)   a(2)
      a(3)     0   -a(1)
     -a(2)   a(1)     0 ];
 