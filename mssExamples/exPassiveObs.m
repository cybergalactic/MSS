% exPassiveObs is compatible with MATLAB and GNU Octave.
% The script plots the loop transfer function of the passive observer in
% for dynamic positioning.

clear;
close all;
clc;

% Filter parameters
lambda = 0.1;
wo     = 0.8976;
wc     = 1.1;
T      = diag([1000 1000 1000]);

% Filter gains
K1 = [ ...
    -2*(1-lambda)*(wc/wo)*eye(3);
     2*wo*(1-lambda)*eye(3) ]

K2 = diag([wc wc wc])
K4 = diag([0.1 0.1 0.01])
K3 = 0.1 * K4

% Mode names
modeName = {'Surge', 'Sway', 'Yaw'};

% Plot all three modes
% Frequency vector
w = logspace(-4,1.5,100);

% Wave response spectrum
wave = tf([wo^2 0],[1 2*lambda*wo wo^2]);

% Store transfer functions
H = cell(3,1);

for i = 1:3

    h0 = tf([1 2*lambda*wo wo^2], ...
        [1 ...
        K1(i+3,i)+K2(i,i)+2*lambda*wo ...
        wo^2+2*lambda*wo*K2(i,i)-K1(i,i)*wo^2 ...
        wo^2*K2(i,i)]);

    hB = tf(K4(i,i)*[1 (1/T(i,i)+K3(i,i)/K4(i,i))], ...
            [1 1/T(i,i)]);

    H{i} = series(h0,hB);

end

figure
bode(H{1},H{2},H{3},wave,w,'LineWidth',2)
grid on

legend('Surge','Sway','Yaw','Wave spectrum','Location','SouthWest','FontSize',12)
