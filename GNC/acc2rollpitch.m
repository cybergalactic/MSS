function [phi, theta] = acc2rollpitch(f)
% [phi, theta] = acc2rollpitch(f) computes the static roll-pitch angles
% phi and theta from 3-axis specific force measurements f = [fx, fy, fz]
%
% Author:    Thor I. Fossen
% Date:      21 March 2020
% Revisions:

phi = atan( f(2) / f(3) );
theta = atan( f(1)/ sqrt( f(2)^2 + f(3)^2 ) );

end