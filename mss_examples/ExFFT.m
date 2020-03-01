% ExFFT Estimation of wave encounter frequency by FFT. Time-domain data are
% generated using a wave spectrum transfer function The FFT is used to
% compute the single-sided spectra for both signal and the peak frequency
% of the spectra are estimated.
%
% Author:    Thor I. Fossen
% Date:      1 March, 2020
% Revisions: 

we = 0.8;                     % peak frequency [rad/s]

fs = 500;                     % Sampling frequency
h = 1/fs;                     % Sampling time
N = 30*60*fs;                 % 30 minutes data
t = (0:N-1)*h;                % Time vector

% Wave spectrum data and sinusoidal (regular) wave data
Kw = 10; lambda = 0.1;                          % wave spectrum data
sys = tf([Kw 0],[1 2*lambda*we we*we]);          
[mag,phase,wout] = bode(sys,logspace(-1,0.3,1000));
mag = reshape(mag(1,:),1,1000);
x1 = lsim(sys,randn(1,length(t)),t,0,'zoh')';   % time responses
x2 = cos(we * t);  
X = [x1; x2];

% Fast Fourier transform (FFT)
n = 2^nextpow2(N);      % Pad the input with trailing zeros
Y = fft(X,n,2);         % Compute the FFT
P2 = abs(Y/N);          % Double-sided spectrum of each signal
P1 = P2(:,1:n/2+1);     % Sigle-sided spectrum of each signal
P1(:,2:end-1) = 2*P1(:,2:end-1);

% Plots
f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors
M = 600;                               % no of samples to plot
figure(gcf);
subplot(2,1,1); 
plot(w(1:M),P1(1,1:M)/max(P1(1,1:M)),wout,mag/max(mag),'linewidth',2);
title(['Normalized wave spectrum in the frequency domain']); grid;
subplot(2,1,2);
plot(w(1:M),P1(2,1:M),'linewidth',2);    
title(['Normalized sinusoidal in the frequency domain']); grid;
