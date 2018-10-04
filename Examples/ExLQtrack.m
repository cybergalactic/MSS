% ExLQtrack  Computes the LQ optimal tracking gains for a mass-damper system
% Author:    Thor I. Fossen
% Date:      25 June 2001
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

% Design matrices
Q = diag([1]);      % user editable tracking error weights
R = diag([1]);      % user editable input weights
%
% System matrices
A = [0 1; -1 -2];   % user editable state matrix
B = [0; 1];         % user editable input matrix
C = [1 0];          % user editable output matrix

% Compute optimal gain matrices
[G1,G2] = lqtracker(A,B,C,Q,R)
