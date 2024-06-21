% exLinspec requires the MATLAB Optimization Toolbox.
% Linear approximation to PM, JONSWAP og Torsethaugen spectra using 
% nonlinear least-squares(NLS).
%
% Author:    Thor I. Fossen
% Date:      2001-08-15 
% Revisions: 
%   2008-12-10  Use updated function wavespec.m
%   2021-07-03  Function Slin.m added at the end of the file 

global sigma wo

wo = 0.8;  
To = 2 * pi / wo;
Hs = 10;
wmax = 3;
w = (0.0001:0.01:wmax)';  

% Modified PM
figure(gcf)
subplot(311)
S = wavespec(3,[Hs,To],w,1);  
sigma = sqrt(max(S));

lambda = lsqcurvefit(@Slin, 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('Modified PM spectrum','Linear approximation')

% JONSWAP
subplot(312)
S = wavespec(7,[Hs,wo,3.3],w,1);   sigma = sqrt(max(S));

lambda = lsqcurvefit(@Slin, 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('JONSWAP spectrum','Linear approximation')

% Torsethaugen (only one peak is fitted)
subplot(313)
S = wavespec(8,[Hs,wo],w,1);   sigma = sqrt(max(S));

lambda = lsqcurvefit(@Slin, 0.1, w, S)

hold on; plot(w,Slin(lambda,w),'linewidth',2); hold off; 
legend('Torsethaugen spectrum','Linear approximation')

set(gca,'FontSize',11)

%% Function Slin
function Pyy = Slin(lambda,w)
% Pyy = Slin(lambda,w) 2nd-order linear power spectracl density (PSD) function
%
% w       : Wave spectrum frequency (rad/s)
% lambda  : Relative damping factor
%
% Author:   Thor I. Fossen
% Date:     2001-08-15 
% Revisions:

global sigma wo

Pyy = 4*(lambda*wo*sigma)^2*w.^2 ./ ( (wo^2-w.^2).^2 + 4*(lambda*wo.*w).^2 );

end

