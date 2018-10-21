function [SVD]=imp2svd(Y)
% Y - is verctor containg the samples of the impulse
%       response of a SISO system.
%
% Author: Tristan Perez 
% Date:   2007.8.4
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
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


[mm,nn] = size(Y);
if nn==1,
  Y=Y';
  [mm,nn] = size(Y);
end
nt=nn-1;
h=hankel(Y(2:nn));
[u,SVDM,v] = svd(h);
SVD=SVDM*ones(max(size(SVDM)),1);