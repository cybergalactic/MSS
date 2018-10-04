function C = m2c(M,nu)
% C = M2C(M,nu) computes the 6x6 Coriolis-centripetal matrix C from the
% the 6x6 system inertia matrix M>0 and the 6x1 velocity vector nu 
% 
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 26th June 2002,  M21 = M12 is corrected to M12'
%            10th Jan 2004, the computation of C is generalized to a nonsymmetric M > 0
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


nu1 = nu(1:3);
nu2 = nu(4:6);

% T = 1/2 (v'Mv) = 1/2 v'(1/2(M + M')+1/2(M-M')) v =  1/2 v' Msym v
% Msym = 1/2 (M+M')

Msym = 0.5*(M+M');
M = Msym;

M11 = M(1:3,1:3);
M12 = M(1:3,4:6);
M21 = M12';
M22 = M(4:6,4:6);

dt_dnu1 = M11*nu1 + M12*nu2;
dt_dnu2 = M21*nu1 + M22*nu2;

C = [  zeros(3,3)      -Smtrx(dt_dnu1)
      -Smtrx(dt_dnu1)  -Smtrx(dt_dnu2) ];
