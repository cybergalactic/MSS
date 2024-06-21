% exPullout is compatible with MATLAB and GNU Octave (www.octave.org).  
% The script performs a pullout maneuver for two different ships.
%
% Author:   Thor I. Fossen
% Date:     25th July 2001
% Revisions: 

delta_c = deg2rad(20);   % Rudder angle for manuver (rad)
h = 0.1;                 % Sampling time (sec)

disp('Pullout maneuver for the Mariner-class cargo vessel (stable ship)')

%% Mariner-class cargo ship, cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);        % x  = [ u v r x y psi delta ]' (initial values)
ui = delta_c;           % ui = delta_c 
 
pullout('mariner',x,ui,h);

%% The Esso Osaka tanker (see tanker.m)
disp(' '); 
disp('Press a key to simulate: The Esso Osaka tanker (unstable ship)'); 
pause

n = 80;
U = 8.23;
x = [ U 0 0 0 0 0 0 n ]';  % x = [ u v r x y psi delta n ]' (initial values)
n_c     = 80;              % n_c = propeller revolution in rpm
depth   = 200;             % Water depth
ui = [delta_c, n_c, depth];
[t,r1,r2] = pullout('tanker',x,ui,h);
