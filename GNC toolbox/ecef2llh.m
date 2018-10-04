function [l,mu,h] = ecef2llh(x,y,z)
% [l,mu,h] = ECEF2LLH(x,y,z) computes the longitude l (rad),
% latitude mu (rad) and height h (m) from the ECEF positions (x,y,z)
%
% Author:   Thor I. Fossen
% Date:     7th June 2001
% Revisions: 1st September 2002, atan2(y/x) replaced by atan2(y,x) 
%            2nd September 2002, new output argument for height h was added
%            27th January, 2003, angle outputs are defined in rad
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>


r_e = 6378137;          % WGS-84 data
r_p = 6356752;
e = 0.08181979099211;
l = atan2(y,x);
eps = 1;
tol = 1e-10;
p = sqrt(x^2+y^2);
mu = atan(z/(p*(1-e^2)));

while (eps > tol),
   N = r_e^2/sqrt(r_e^2*cos(mu)^2+r_p^2*sin(mu)^2);
   h = p/cos(mu)-N;
   mu0 = mu;
   mu = atan(z/(p*(1-e^2*N/(N+h))));
   eps = abs(mu-mu0);
end

%mu = mu*180/pi;
%l = l*180/pi;
