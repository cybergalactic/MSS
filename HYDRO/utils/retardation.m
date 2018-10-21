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
