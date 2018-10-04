function  [Vw]=hs2vw(Hs)
% This function converts the significan wave heihgt Hs into an equivalent wind
% speed Vw accodding to 
% 
%        Vw= sqrt((Hs*g/0.21));
%
% Use:  [Vw]=hs2vw(Hs)
% Hs -Significant wave height [m].
% Vw -Wind speed [m/s]
%
% Reference:
% T.I. Fossen (2002) "Marine Control Systems" Marine Cybernetics. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by: Tristan Perez in 2005  
% Last mod. by: TP
% Date: 2005-03-12
% Comments:  created for support of the old function pierson.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


g=9.81;
Vw= sqrt((Hs*g/0.21));