function S = wavespec(spectrumType, Parameter, W, PlotFlag)
% S = wavespec(spectrumType, Parameter, W, PlotFlag) is obsolete and has 
% been replaced by S = waveSpectrum(spectrumType, Parameter, W, PlotFlag).
%
% Author:    Thor I. Fossen
% Date:      2024-06-07

disp('S = wavespec(spectrumType, Parameter, W, PlotFlag) is obsolete.')
disp('Use S = waveSpectrum(spectrumType, Parameter, W, PlotFlag)')

S = waveSpectrum(spectrumType, Parameter, W, PlotFlag);

end