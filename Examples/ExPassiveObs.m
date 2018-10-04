% ExPassiveObs Plots the loop transfer function of the passive observer
% Author:      Thor I. Fossen
% Date:        15th October 2002
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

% filter paramaters
T = eye(3)/1000
lambda = 0.1;
wo   = 0.8976;
wc   = 1.1;
T = diag([1000 1000 1000]);

% filter gains
K1 = [-2*(1-lambda)*(wc/wo)*eye(3)
    2*wo*(1-lambda)*eye(3)     ]
K2 = diag([wc wc wc])
K4 = diag([0.1 0.1 0.01])
K3 = 0.1*K4

% Bode plots
figure(gcf)
w = logspace(-4,1.5,100);

i=1;  % mode = surge
h0 = tf([1 2*lambda*wo wo*wo],[1 (K1(i+3,i)+K2(i,i)+2*lambda*wo) (wo^2+2*lambda*wo*K2(i,i)-K1(i,i)*wo^2) wo^2*K2(i,i)]) 
hB = tf(K4(i,i)*[1 (1/T(i,i)+K3(i,i)/K4(i,i))],[1 1/T(i,i)])
wave = tf([wo^2 0],[1 2*lambda*wo wo^2])      % wave response spectrum

bode(series(h0,hB),wave,w)





