% ExLinspec  linear approx. to PM, JONSWAP og Torsethaugen spectra using 
%            nonlinear least-squares(NLS)
%
% Author:    Thor I. Fossen
% Date:      15.08.2001
% Revisions: 10.12.2008  use updated function wavespec.m

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


