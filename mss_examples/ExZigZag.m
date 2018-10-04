% ExZigZag   Generates the zigzag maneuver for two different ships
%
% Author:   Thor I. Fossen
% Date:     23th July 2001
% Revisions: 

t_final = 600;           % final simulation time (sec)
t_rudderexecute = 10;    % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

disp('20-20 Zigzag maneuver for the Mariner class vessel')

% 20-20 zigzag maneuver for the Mariner class cargo ship
% cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = 0;            % delta_c = 0  for time t < t_rudderexecute
 
[t,u,v,r,x,y,psi,U] = zigzag('mariner',x,ui,t_final,t_rudderexecute,h,[20,20]);

% 20-10 zigzag maneuver for a container ship 
% cruise speed 8.0 m/s see container.m)
disp(' '); disp('Press a key to simulate: 20-10 zigzag manuver for a container ship'); pause
 
x     = [8.0 0 0 0 0 0 0 0 0 70]';     % x = [ u v r x y psi delta n ]' (initial values)
delta_c = 0;                           % delta_c = 0 for time t < t_rudderexecute
n_c     = 80;                          % n_c = propeller revolution in rpm
 
ui = [delta_c, n_c];
[t,u,v,r,x,y,psi,U] = zigzag('container',x,ui,t_final,t_rudderexecute,h,[20,10]);
