% ExRefMod   2nd-order reference model with nonlinear damping and velocity saturation
% Author:    Thor I. Fossen
% Date:      3rd November 2001
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

z = 1;     % relative damping ratio
w = 1;     % natural frequency
delta = 1; % nonlinear damping coeff.
vmax = 1;  % max velocity
h = 0.1;   % sampling time
N = 200;   % number of samples

t  = 0:h:h*N;

% linear mass-damper-spring system
r1 = 10*ones(max(size(t)),1);
[A,B,C,D] = ord2(w,z); 
[x1,y1]   = lsim(A,B,C,D,r1,t); 

r2 = 10*r1;
[A,B,C,D] = ord2(w,z); 
[x2,y2]   = lsim(A,B,C,D,r2,t); 

% nonlinear damping
x = 0;
v = 0;
r2 = 10;
y2 = zeros(N+1,2);
for i=1:N+1,
   y2(i,:) = [x v];   
   x_dot = v;
   v_dot = w^2*(r2-x) - 2*z*w*v - delta*abs(v)*v;
   v = v + h*v_dot;
   x = x + h*x_dot;
end

% velocity saturation
x = 0;
v = 0;
r3 = 10;
y3 = zeros(N+1,2);
for i=1:N+1,
   y3(i,:) = [x v];   
   v_dot = w^2*(r3-x) - 2*z*w*v;
   x_dot = v;
   v = v + h*v_dot;
   if abs(v)>vmax,       % saturation
      v = sign(v)*vmax;
   end
   x = x + h*x_dot;
end

% plots
figure(gcf)
subplot(211); plot(t,y1(:,1),'linewidth',2)
hold on; plot(t,y2(:,1),'--k',t,y3(:,1),'r'); hold off; grid
title('2nd-order mass-damper-spring reference model')
legend('linear damping','nonlinear damping','velocity saturation')
subplot(212); plot(t,y1(:,2),'linewidth',2);
hold on; plot(t,y2(:,2),'--k',t,y3(:,2),'r'); hold off; grid


