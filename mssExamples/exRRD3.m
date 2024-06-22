% exRRD3: Inverse response in roll due to right-half-plane zero 
% (non-minimum phase) for the Son and Nomoto (1982) container ship.
%
% References:  
%   K. Son og K. Nomoto (1982). On the Coupled Motion of Steering and 
%   Rolling of a High Speed Container Ship, Naval Architect of 
%   Ocean Engineering 20:73-83. From J.S.N.A., Japan, Vol. 150, 1981.
%
% Author:    Thor I. Fossen
% Date:      2001-10-31
% Revisions: 

conversion % Load conversion factors

u0     = 7.3;  % Initial speed (m/s)
n_ref =  70;   % Desried propeller speed (RPM)

N = 10000;     % No. of samples
h = 0.05;      % Sample time in seconds

% x = [u v r x y psi p phi delta n ]'
x = [u0 0 0 0 0 0 0 0 0 70]';     
U = u0;

xout = zeros(N+1,length(x)+1);

for i=1:N+1
    
    % Rudder commands
    if (i > 1000 && i < 1500)
        u_ref = deg2rad(10);
    else
        u_ref = 0;
    end
        
    xout(i,:) = [x',U];
    
    % Nonlinear model of Son and Nomoto (1982)
    [xdot,U] = container(x,[u_ref; n_ref]); 
    
     x = x + h*xdot;
end

t = h*(0:1:N)';
figure(gcf)

subplot(311)
plot(t, rad2deg(xout(:,6)), t, 10*rad2deg(xout(:,8)))
title('Roll angle and yaw angle for a 50 s step in rudder angle (deg)')
grid
legend('10\cdot\phi (roll angle)','\psi (yaw angle)')

subplot(312)
plot(t, rad2deg(xout(:,9)))
title('Rudder angle (deg)')
grid

subplot(313)
plot(t, xout(:,11))
title('Ship speed (m/s)')
grid

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
