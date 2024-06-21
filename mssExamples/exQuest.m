% exQuest is compatible with MATLAB and GNU Octave (www.octave.org).
% 6-DOF position/attitude vector from camera measurements using the QUEST 
% (QUaternion ESTimator) algorithm. The QUEST algorithm is a widely used 
% method for attitude determination, particularly for spacecraft. It 
% estimates the orientation of a body with respect to an inertial frame 
% of reference using vector observations and the quaternion representation 
% of rotations. The algorithm is particularly efficient and robust, making 
% it a preferred choice for real-time applications; see quest6dof.m and 
% quest.m.   
%
% References: 
%   Schuster, M. D., Oh, S. D., 1981, Three-axis attitude determination 
%   from vector observations, Journal of Guidance, Dynamics and Control 
%   JGC-4(1):70-77.
%
% Author:    Thor I. Fossen and K. P. Lindegaard
% Date:      5 August 2001
% Revisions: 

% Marker positions
mb1 = [0.69  0  -0.15]';
mb2 = [0.42  0  -0.51]';
mb3 = [-0.38 0  -0.18]';

% Camera position
rcamera = [2.11  -0.48  1.89]';

% Generate camera measurement for ship position: rcg = [ 2.5 3.0 0.1]' and 
% attitude: phi = 10 deg, theta = 5 deg, and psi = 134 deg.
disp('True ship position and attitude:')
rcg = [ 2.5 3.0 0.1]'                        % Location of CG
phi = 10
theta = 5
psi = 134
q = euler2q(phi*pi/180,theta*pi/180,psi*pi/180);  % Quaternions
R = Rquat(q);                                     % Rotation matrix
disp('Camera measurements:')
y1 = rcg - rcamera + R*mb1                        % Camera measurements
y2 = rcg - rcamera + R*mb2
y3 = rcg - rcamera + R*mb3

%% QUEST ALGORITHM
% Compute position/attitude eta in 6DOF, quaternions q and R from the 
% camera measurements y1,y2,y3
y  = [y1; y2; y3];
mb = [mb1; mb2; mb3];
[eta,q,R] = quest6dof(y,mb,rcamera);
[phi,theta,psi] = q2euler(q);

disp('6-DOF position computed from camera measurements using the QUEST algorithm:')
rcg_quest    = eta(1:3)'
phi_quest    = rad2deg(phi)
theta_quest  = rad2deg(theta)
psi_quest    = rad2deg(psi)
