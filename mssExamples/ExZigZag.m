% ExZigZag   Zigzag maneuvers for the Mariner class vessel, a container 
%            ship and the Remus 100 AUV.
%
% Author:    Thor I. Fossen
% Date:      23 Jul 2001
% Revisions: 16 Nov 2023 - included the Remus 100 AUV

t_final = 600;           % final simulation time (sec)
t_rudderexecute = 10;    % time rudder is executed (sec)

disp('20-20 zigzag maneuver for the Mariner class vessel')

% 20-20 zigzag maneuver for the Mariner class cargo ship (mariner.m)
% cruise speed U0 = 7.7 m/s 
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = 0;            % delta_c = 0  for time t < t_rudderexecute

zigzag('mariner',x,ui,t_final,t_rudderexecute,0.1,[20,20]);

% 20-10 zigzag maneuver for the  container ship (container.m)
% Cruise speed 8.0 m/s
disp(' '); 
disp('Press a key to simulate: 20-10 zigzag manuver for a container ship'); 
pause;
 
x     = [8.0 0 0 0 0 0 0 0 0 70]';   % x = [ u v r x y psi delta n ]' (initial values)
delta_c = 0;                         % delta_c = 0 for time t < t_rudderexecute
n_c     = 80;                        % n_c = propeller speed in rpm
ui = [delta_c, n_c];

zigzag('container',x,ui,t_final,t_rudderexecute,0.1,[20,10]);

% 10-10 zigzag maneuver for the Remus 100 AUV (remus100.m)
disp(' '); 
disp('Press a key to simulate: 20-20 zigzag manuver for the Remus100 AUV'); 
pause;

t_final = 100;           % final simulation time (sec) 
x       = zeros(12,1);   % initial state vector
delta_r = 0;             % delta_r = 0 for time t < t_rudderexecute
delta_s = 0;             % stern dive plane
n_c     = 300;           % n_c = propeller speed in rpm 
ui = [delta_r, delta_s, n_c];

zigzag6dof('remus100',zeros(12,1),ui,t_final,t_rudderexecute,0.02,[20,20]);


