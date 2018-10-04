function [We]=ww2we(chi,U,Ww)
% Function to transform from a vector wave frequancy to 
% encounter frequencies
%
%Use: [We]=ww2we(chi,U,W)
% 
% We   vecrtor of encounter fequency values [rad/sec] 
% chi  encounter angle [rad] 0-following seas, pi-head seas
% U    forward speed [m/sec]
% W    vecrtor of wave fequancy values [rad/sec] 
%
%Reference: A.R.M.J LLoyd "Seakeeping: Ship Behaviour in Rough Wheather."
%John Wiley & Sons, 1989.
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

w=[];
for k=1:length(Ww)
      we_aux=(Ww(k)-Ww(k)^2*U*cos(chi)/9.81);
   w=[w;we_aux];
end
We=w;