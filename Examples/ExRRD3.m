% ExRRD3     Inverse response in roll due to right-half-plane zero (non-minimum phase)
%            Son og Nomoto container ship
% Author:    Thor I. Fossen
% Date:      31 October 2001
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

conversion % load conversion factors

u0     = 7.3;    % initial speed
n_ref =  70;     % rpm

N = 10000;  % no. of samples
h = 0.05;   % sample time

x = [u0 0 0 0 0 0 0 0 0 70]';     % x = [u v r x y psi p phi delta n ]'
U = u0;

xout = zeros(N+1,length(x)+1);

for i=1:N+1,
    
    % 10 deg rudder command
    if (i>1000 & i<1500),
        u_ref = 10*D2R;
    else
        u_ref = 0;
    end
        
    xout(i,:) = [x',U];
    
    % nonlinear model of Son and Nomoto
    [xdot,U] = container(x,[u_ref; n_ref]); 
    
     x = x + h*xdot;
end

t = h*(0:1:N)';
figure(gcf)

subplot(211)
plot(t,R2D*xout(:,6),t,10*R2D*xout(:,8))
title('Roll angle 10\cdot\phi (deg) and yaw angle \psi (deg) for a 10 (deg) rudder step in 50 (s)')
grid

subplot(212)
plot(t,xout(:,11))
title('Forward speed (m/s)')
grid