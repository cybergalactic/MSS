 echo off
% KINDEMO   Demonstration of kinematic computations. Comparative study of
% Euler angles and unit quaternions.
%
% Author:    Thor I. Fossen
% Date:      14 Jun 2014
% Revisions: 07 May 2022 - changed to exact discretization of quaternions

echo on; clc

h = 0.05;     % sampling time
N = 200;     % number of samples

D2R = pi/180;   
R2D = 180/pi;

phi   = 10 * D2R;   % initial Euler angles
theta =  5 * D2R; 
psi   =  1 * D2R;

nu = sin((0:N)*h)'*[0.5 0.1 -0.2];  % angular velocity vector

pause % Strike any key to continue

%##########################################################################
% Euler angle kinematics
%##########################################################################
x = [phi theta psi]';

xdata = zeros(N,7);
for i = 1:N
   q = euler2q(x(1),x(2),x(3)); % transform the Euler angles to a quaternion
   xdata(i,:) = [q' x'];        % store data
   
   [J,J1,J2] = eulerang(x(1),x(2),x(3));  % transformation matrices 
   dx = J2*nu(i,:)';             % Euler angle differential equation
   x  = euler2(dx,x,h);          % Euler's integration method
   echo off
end
echo on; 

%##########################################################################
% Quaternion kinematics
%##########################################################################
q = euler2q(phi,theta,psi); 

qdata = zeros(N,7);
for i = 1:N
   [phi,theta,psi] = q2euler(q);     % transform q to Euler angles
   qdata(i,:) = [q' phi theta psi];  % store data
   q = expm(Tquat(nu(i,:)')*h)*q;    % exact discretization
   
   % Euler's integration method (alternative to exact discretization)
   % [J,J1,J2] = quatern(q);        % transformation matrices
   % dq = J2*nu(i,:)';              % quaternion differential equation
   % q  = euler2(dq,q,h);              
  
   q  = q/norm(q);                  % normalization of q 
   echo off
end
echo on;  

pause % Strike any key to continue
echo off

%##########################################################################
% Plots
%##########################################################################
t = (0:h:(N-1)*h)';
figure(1)
subplot(221); plot(t,qdata(:,1:4),'linewidth',2); 
title('Unit quaternion q from ODE','fontsize',14); 
xlabel('time (s)');  grid;
subplot(222); plot(t,R2D*qdata(:,5:7),'linewidth',2); 
title('[\phi,\theta,\psi] = q2euler(q)','fontsize',14); 
xlabel('time (s)');  grid;
subplot(223); plot(t,xdata(:,1:4),'linewidth',2); 
title('q = euler2q(\phi,\theta,\psi)','fontsize',14); 
xlabel('time (s)'); ylabel('deg'); grid;
subplot(224); plot(t,R2D*xdata(:,5:7),'linewidth',2);
title('Euler angles \phi,\theta and \psi from ODE','fontsize',14); 
xlabel('time (s)'); ylabel('deg'); grid;

%##########################################################################
echo on
pause % Strike any key for main menu
echo off
disp('End of demo')

