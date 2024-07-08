% exLinspec requires the MATLAB Optimization Toolbox.
% Linear approximation to PM, JONSWAP, and Torsethaugen spectra using 
% nonlinear least-squares (NLS) curve fitting.
%
% Author:    Thor I. Fossen
% Date:      2001-08-15 
% Revisions: 
%   2008-12-10  Use updated function wavespec.m
%   2021-07-03  Function Slin.m added at the end of the file 
%   2024-07-07  Use new function waveSpectrum.m. Removed the global variable

clearvars;

wo = 0.8;  
To = 2 * pi / wo;
Hs = 10;
wmax = 3;
w = (0:0.01:wmax)';  

% Initial guess for lambda
lambda0 = 0.1;

figure(gcf); clf;

% Modified PM
subplot(311)
S = waveSpectrum(3, [Hs, To], w, 1);  
sigma = sqrt(max(S));
lambda = lsqcurvefit(@(lambda, w) Slin(lambda, wo, sigma, w), lambda0, w, S);
disp(['lambda_PM = ', num2str(lambda)])
hold on; plot(w, Slin(lambda, wo, sigma, w), 'b'); hold off; grid;
legend('Modified PM spectrum', 'Linear approximation')

% JONSWAP
subplot(312)
S = waveSpectrum(7, [Hs, wo, 3.3], w, 1);   
sigma = sqrt(max(S));
lambda = lsqcurvefit(@(lambda, w) Slin(lambda, wo, sigma, w), lambda0, w, S);
disp(['lambda_JONSWAP = ', num2str(lambda)])
hold on; plot(w, Slin(lambda, wo, sigma, w), 'b'); hold off; grid;
legend('JONSWAP spectrum', 'Linear approximation')

% Torsethaugen (only one peak is fitted)
subplot(313)
S = waveSpectrum(8, [Hs, wo], w, 1);   
sigma = sqrt(max(S));
lambda = lsqcurvefit(@(lambda, w) Slin(lambda, wo, sigma, w), lambda0, w, S);
disp(['lambda_Torsethaugen = ', num2str(lambda)])
hold on; plot(w, Slin(lambda, wo, sigma, w), 'b'); hold off; grid; 
legend('Torsethaugen spectrum', 'Linear approximation')

set(findall(gcf, 'type', 'line'), 'linewidth', 2)
set(findall(gcf, 'type', 'text'), 'FontSize', 12)
set(findall(gcf, 'type', 'legend'), 'FontSize', 12);

%% Function Slin
function Pyy = Slin(lambda, wo, sigma, w)
% Pyy = Slin(lambda, wo, sigma, w) 2nd-order linear power spectral density 
% (PSD) function.
%
% lambda  : Relative damping factor
% wo      : Peak frequency (rad/s)
% sigma   : Spectrum gain
% w       : Wave spectrum frequency (rad/s)
%
% Author:   Thor I. Fossen
% Date:     2001-08-15 
% Revisions:

Pyy = 4 * (lambda * wo * sigma)^2 * w.^2 ./ ...
    ( (wo^2 - w.^2).^2 + 4 * (lambda * wo .* w).^2 );

end