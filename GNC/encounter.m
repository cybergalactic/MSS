function We = encounter(Ww,U,beta)
% We = ENCOUNTER(Ww,U,beta) computes the encounter frequency We (rad/s) as a function
%       of wave peak frequeny Ww (rad/s), vessel speed U (m/s) and wave direction 
%       beta, 0 (deg) for following seas and 180 (deg) for head seas.
%       
% We = abs(Ww - Ww.^2*U*cos(beta*pi/180)/g);
% size(We)=size(Ww)
%
% Author:   Thor I. Fossen
% Date:     2nd November 2001
% Revisions: 

g = 9.81; % acceleration of gravity
We = abs(Ww - Ww.^2*U*cos(beta*pi/180)/g);