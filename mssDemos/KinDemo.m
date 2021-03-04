echo off
% KINDEMO   Demonstration of kinematic computations
%
% Author:   Thor I. Fossen
% Date:     2001-06-14
% Revisions: 

echo on; clc

h = 0.1;    % samling time
N = 200;    % number of samples

deg2rad = pi/180;   rad2deg = 180/pi;

phi   = 10*deg2rad; % initial Euler angles
theta =  5*deg2rad; 
psi   =  1*deg2rad;

nu = sin((0:N)*h)'*[0.5 0.1 -0.2];  % angular velocity vector

pause % Strike any key to continue

%#################################################################
% Euler angle kinematics
%#################################################################
x = [phi theta psi]';

xdata = zeros(N+1,7);
for i = 1:N
   q = euler2q(x(1),x(2),x(3));  % transform Euler angles to q
   xdata(i,:) = [q' x'];         % store data
   
   [J,J1,J2] = eulerang(x(1),x(2),x(3));  % transformation matrices 
   dx = J2*nu(i,:)';             % Euler angle diff. eq.
   x  = euler2(dx,x,h);          % Euler integration
   echo off
end
echo on; 

%#################################################################
% Quaternion kinematics
%#################################################################
q = euler2q(phi,theta,psi); 

qdata = zeros(N+1,7);
for i = 1:N
   [phi,theta,psi] = q2euler(q);     % transform q to Euler angles
   qdata(i,:) = [q' phi theta psi];  % store data
   
   [J,J1,J2] = quatern(q);  % transformation matrices
   dq = J2*nu(i,:)';        % quaternion diff. eq.
   q  = euler2(dq,q,h);     % Euler integration
   q  = q/norm(q);          % normalization of q 
   echo off
end
echo on;  

pause % Strike any key to continue
echo off

%#################################################################
% Plots
%#################################################################
t = (0:h:N*h)';
figure(1)
subplot(221); plot(t,qdata(:,1:4)); 
title('q from ODE'); xlabel('time (s)'); 
subplot(222); plot(t,rad2deg*qdata(:,5:7)); 
title('[\phi,\theta,\psi] = q2euler(q)'); xlabel('time (s)'); 
subplot(223); plot(t,xdata(:,1:4)); 
title('q = euler2q(\phi,\theta,\psi)'); xlabel('time (s)'); ylabel('deg')
subplot(224); plot(t,rad2deg*xdata(:,5:7)); 
title('\phi,\theta,\psi from ODE'); xlabel('time (s)'); ylabel('deg')
%#################################################################

echo on
pause % Strike any key for main menu
close(1); echo off
disp('End of demo')

