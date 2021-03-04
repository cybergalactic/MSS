% ExTurnCircle  Generates the turning circle for two different ships
%
% Author:   Thor I. Fossen
% Date:     21th July 2001
% Revisions: 11th June 2003 (Thor I. Fossen) rudder execute is changed 
%                                 from 15 deg to 35 deg (IMO regulations)

t_final = 700;          % final simulation time (sec)
t_rudderexecute = 100;   % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

disp('Turning cirlce for the Mariner class vessel')

% Mariner class cargo ship, cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = -35*pi/180;   % delta_c = -delta_R at time t = t_rudderexecute
 
[t,u,v,r,x,y,psi,U] = turncircle('mariner',x,ui,t_final,t_rudderexecute,h);

% Container ship (see container.m)
disp(' '); disp('Press a key to simulate: Container Ship'); pause

x     = [8.0 0 0 0 0 0 0 0 0 70]';     % x = [ u v r x y psi delta n ]' (initial values)
delta_c = -35*pi/180;                  % delta_c = -delta_R at time t = t_rudderexecute
n_c     = 80;                          % n_c = propeller revolution in rpm

ui = [delta_c, n_c];
[t,u,v,r,x,y,psi,U] = turncircle('container',x,ui,t_final,t_rudderexecute,h);
