function [ph]=rand_phases(nc)
%Function to generate a uniformely 
%distributed vector of random phases 
%in the interval [-pi pi].
%
% Use: [ph]=rand_phases(nc)
%
% ph -vector of random phases 
% nc- number of components
% 
% Created by: Tristan Perez in 2001  
% Last mod. by: Tristan Perez 
% Date: 9 March 2005
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



ph=2*pi*(rand([nc 1])-.5);
