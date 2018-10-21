function Kij = retardation(t,w,B,i,j)
% Kij = retardation(t,w,B,i,j) computes the retardation function:
%
%   Kij(t) = int_0^t Bij(w)*cos(w*tau)*dtau
%
% using the trapezoidal rule (w must be equally spaced). For SISO systems
% use K = retardation(t,w,B)
%
% Thor I. Fossen
% 2005-05-29
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

K = zeros(1,length(t));
dw = w(2)-w(1);

if nargin ==3,

    for n1=1:length(t),
        K(n1) = 2/pi/2*dw*( B(1) + B(length(w))*cos(w(length(w))*t(n1)) );
        for n2=2:(length(w)-1),
            K(n1) = K(n1) + 2/pi*dw*B(n2)*cos(w(n2)*t(n1));
        end
    end

elseif nargin == 5,

    for n1=1:length(t),
        K(n1) = 2/pi/2*dw*( B(i,j,1) + B(i,j,length(w))*cos(w(length(w))*t(n1)) );
        for n2=2:(length(w)-1),
            K(n1) = K(n1) + 2/pi*dw*B(i,j,n2)*cos(w(n2)*t(n1));
        end
    end

end

Kij = K;
