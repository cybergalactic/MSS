function [phi, theta] = acc2rollpitch(f)
% acc2rollpitch is compatible with MATLAB and GNU Octave (www.octave.org). 
% The function [phi, theta] = acc2rollpitch(f) computes the static 
% roll-pitch angles phi and theta from 3-axis specific force measurements 
% f = [fx, fy, fz].
%
% Author:    Thor I. Fossen
% Date:      2020-03-20
% Revisions:

phi = atan2( f(2), f(3) );  
theta = atan2( f(1), sqrt( f(2)^2 + f(3)^2 ) ); 

end
