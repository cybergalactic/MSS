function H = Hmtrx(r)
% H = HMTRX(r) computes the 6x6 system transformation matrix
% 
% S = Smtrx(r);
% H = [eye(3)     S'
%      zeros(3,3) eye(3) ];       Property: inv(H(r)) = H(-r)
%
% If r = r_g is the vector from CO to CG, the model matrices in CO and CG
% are related by: M_CO = H(r_g)'*M_CG*H(r_g). Generalized position and
% force satisfy: eta_CO = H(r_g)'*eta_CG and tau_CO = H(r_g)'*tau_CG 
%
% Author:    Thor I. Fossen
% Date:      2001-06-14
% Revisions: 2008-01-17  R1.0
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


S = Smtrx(r);
H = [eye(3)     S'
     zeros(3,3) eye(3) ];
