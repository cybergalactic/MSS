function [G1,G2,G3] = lqtracker(A,B,C,Q,R,E)
% LQTRACKER  computes the LQ tracker gain matrices for LTI systems:
%     dx/dt = Ax + Bu + Ew where E is an optionally input for diturbance feedforward
%     [G1,G2]    = lqtracker(A,B,C,Q,R)   returns u = G1*x + G2*yd 
%     [G1,G2,G3] = lqtracker(A,B,C,Q,R,E) returns u = G1*x + G2*yd + G3*w
%
% Author: Thor I. Fossen
% Date: 16th June 2001
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


[K,P,EE] = lqr(A,B,C'*Q*C,R);
G1 = -inv(R)*B'*P;
Temp = inv((A+B*G1)');
G2 = -inv(R)*B'*Temp*C'*Q;

if nargin==6,
   G3 = inv(R)*B'*Temp*P*E;
else
   G3 = NaN;
end
