% ExLinspec  linear approx. to PM, JONSWAP og Torsethaugen spectra using 
%            nonlinear least-squares(NLS)
%
% Author:    Thor I. Fossen
% Date:      15.08.2001
% Revisions: 10.12.2008  use updtaed function wavespec.m
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org
%
global sigma wo

wo = 0.8;  To = 2*pi/wo;
Hs = 10;
wmax = 3;
w = (0.0001:0.01:wmax)';  

% Modified PM
subplot(311)
S = wavespec(3,[Hs,To],w,1);  sigma = sqrt(max(S));

lambda = lsqcurvefit('Slin', 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('Modified PM spectrum','Linear approximation')

% JONSWAP
subplot(312)
S = wavespec(7,[Hs,wo,3.3],w,1);   sigma = sqrt(max(S));

lambda = lsqcurvefit('Slin', 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('JONSWAP spectrum','Linear approximation')

% Torsethaugen (only one peak is fitted)
subplot(313)
S = wavespec(8,[Hs,wo],w,1);   sigma = sqrt(max(S));

lambda = lsqcurvefit('Slin', 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('Torsethaugen spectrum','Linear approximation')

set(gca,'FontSize',11)


