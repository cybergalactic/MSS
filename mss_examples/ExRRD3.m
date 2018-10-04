% ExRRD3     Inverse response in roll due to right-half-plane zero (non-minimum phase)
%            Son og Nomoto container ship
% Author:    Thor I. Fossen
% Date:      31 October 2001
% Revisions: 

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