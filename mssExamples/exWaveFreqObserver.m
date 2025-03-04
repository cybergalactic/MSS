% This script demonstrates the implementation of a nonlinear wave frequency 
% observer for estimating the wave encounter frequency from an INS-measured 
% position or orientation signal. The method is based on the nonlinear observer:
%
% Reference:
%   D. J. Belleter, R. Galeazzi and T. I. Fossen (2015). Experimenal Verification
%   of a Globally Exponentially Stable Nonlinear Wave Encounter Frequency
%   Estimator. Ocean Engineering 97(15), 48–56.
%
% The observer consists of:
%   1. A high-pass filter to remove low-frequency components.
%   2. A nonlinear observer to estimate the wave encounter frequency.
%
% The test signal is a combination of a low-frequency sawtooth wave 
% and a wave-frequency component simulated with white noise.
%
% Author: Thor I. Fossen
% Date: 2025-03-03
% Revisions: 

clearvars; clf; 

% Sampling parameters
h = 0.05; % Sampling time 
T_final = 200; % Final simulation time in seconds

% Time vector initialization
t = (0:h:T_final)'; % Time vector from 0 to T_final          
nTimeSteps = length(t); % Number of time steps

% Wave frequency observer (Belleter, Galeazzi and Fossen, 2015)
K_f = 1; % Observer gain, typically K_max * A_max is 1 to 50
w_f = 4; % Filter natural frequency, larger than the encounter frequency w_e 
x = [0; 0; 0.5^2]; % Initial state vector for the observer
x_hp = 0; % Initial state vector for HP filter

% Test signal
omega_0 = 1.0; % Wave spectrum peak frequency
lambda = 0.1; % Wave spectrum relative damping ratio
K_wave = 1; % Wave amplitude gain

white_noise = randn(size(t)); 
s = tf('s'); 

H_s = K_wave * s / (s^2 + 2 * omega_0 * lambda * s + omega_0^2);
H_z = c2d(H_s, h, 'tutsin');
u_LF = 1 * sawToothWave(0.1 * t);
u = u_LF + lsim(H_z, white_noise, t); % LF + WF position

% Nonlinear wave frequency observer
omega_e = zeros(nTimeSteps,1);
for k = 1:nTimeSteps
    [x_hp, y_obs] = highPassFilter(x_hp, u(k), 0.3, h);
    [x, omega_obs] = waveFreqObserver(x, y_obs, w_f, K_f, h);
    omega_e(k) = omega_obs;
end

%% Plot Results
figure(1); figure(gcf)
subplot(211)
yline(omega_0,'k','linewidth',2); hold on;
plot(t, omega_e, 'b');
hold off
legend('True encounter frequency', 'Estimated Wave Frequency');
xlabel('Time [s]'); 
grid on;
title('Wave Frequency Estimation using Nonlinear Observer');
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(t, u, 'g'); hold on;
plot(t, u_LF, 'k'); 
hold off
legend('LF+WF Test Signal', 'Sawtooth LF Signal');
xlabel('Time [s]'); 
grid on;
title('Position measurement');
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
