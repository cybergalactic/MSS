function [x,y,z] = llh2ecef(l,mu,h)
% [x,y,z] = LLH2ECEF(l,mu,h) computes the  ECEF positions (x,y,z)
% from longitude l (rad), latitude mu (rad) and height h
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 27th January 2003, inputs l and mu are defined in rad
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


r_e = 6378137; % WGS-84 data
r_p = 6356752;
e = 0.08181979099211;
N = r_e^2/sqrt( (r_e*cos(mu))^2 + (r_p*sin(mu))^2 );
x = (N+h)*cos(mu)*cos(l);
y = (N+h)*cos(mu)*sin(l);
z = (N*(r_p/r_e)^2 + h)*sin(mu);