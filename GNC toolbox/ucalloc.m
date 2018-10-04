function u = ucalloc(K,T,W,tau)
% u = ucalloc(K,T,W,tau) unconstrained control allocation. The generalized
%     force vector tau = T*K*u (dim n) is distributed to the input vector u
%     (dim r) where r>=n by minimizing the force f=K*u.
%
%     An unconstrained solution u = inv(K)*inv(W)*T'*inv(T*inv(W)*T')
%     exists if T*T' is non-singular. 
%
%     - K is a diagonal rxr matrix of force coeffisients 
%     - T is a nxr constant configuration matrix. 
%     - W is a rxr positive diagonal matrix weighting (prizing) the 
%       different control forces f = K*u.
%
% Author:    Thor I. Fossen
% Date:      3rd November 2001
% Revisions: 25th September 2002 - function name changed from alloc.m to ucalloc.m
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

if det(T*T')==0, 
    error('T*T''is singular'); 
elseif det(W)==0,
    error('W must be positive'); 
else
    Winv = inv(W);
    u = inv(K)*Winv*T'*inv(T*Winv*T')*tau;
end


