% ExObsCtr   Observability and Controllability of Ships (see supply.m)
% Author:    Thor I. Fossen
% Date:      11 July 2002
% Revisions: 20 Dec 2009   minor modifications
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

format compact
disp('Controllabilty and observability of offShore supply vessel length 76 m')

% Normalization variables

L   =  76.2;           % length of ship (m)
g   =  9.8;            % acceleration of gravity (m/s^2)
m   = 6000e3;          % mass (kg)

T    = diag([1 1 L]);
Tinv = diag([1 1 1/L]);

% Model matricses
Mbis = [1.1274         0          0
             0    1.8902    -0.0744
             0   -0.0744     0.1278];

Dbis = [0.0358        0        0
             0        0.1183  -0.0124
             0       -0.0041   0.0308];
   
Athr = diag([-1/100 -1/100 -1/100]);

M = m*Tinv^2*(T*Mbis*Tinv)
D = m*sqrt(g/L)*Tinv^2*(T*Dbis*Tinv)

% state space model
A = [ zeros(3,3) eye(3)     zeros(3,3)
      zeros(3,3) -inv(M)*D  inv(M)
      zeros(3,3) zeros(3,3) Athr      ];

B = [zeros(3,3); zeros(3,3); -Athr ];

F = [ zeros(3,3) zeros(3,3) eye(3)
      zeros(3,3) zeros(3,3) zeros(3,3)  
      zeros(3,3) inv(M) -inv(M)*D ];

H = [ eye(3) zeros(3,3) zeros(3,3)];


n=rank(ctrb(A,B))

n=rank(obsv(F,H))

% augmented state space modelwith integral action
C = [ eye(3) zeros(3,3) zeros(3,3)];
Aa = [zeros(3,3) C
      zeros(9,3) A ];

Ba = [zeros(3,3)
      B ];

n=rank(ctrb(Aa,Ba))
