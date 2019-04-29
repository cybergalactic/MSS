% ExResonance computes the responses in heave, roll and pitch for a marine
% craft exposed to a regular wave. The closed-form solution of a linear
% mass-damper-spring system with sinusoidal forcing F sin(wt) is used. 
% Harmonic oscillator:
%   ..     .
% m x  + d x + k x = F cos(wt),     w_n = sqrt(k/m)
%
% The results are plotted as change in amplitudes with respect to the 
% frequency ration w/w_n.
%
% Author:    Thor I. Fossen
% Date:      2019-04-29
% Revisions: 

load supply;		% load supply data (ShipX)

% gain for wave-induced force
K = 0.1;    % Force = K * mass * w_n^

% compute relative damping rations and natural frequencies
[T,zeta] = DPperiods(vessel);
z3 = zeta(3);   z4 = zeta(4);   z5 = zeta(5);
w3 = 2*pi/T(3); w4 = 2*pi/T(4); w5 = 2*pi/T(5);

wr = 0:0.01:2.5;   % wr = w / w_n

figure(gcf)

w = w3 * wr;
A3 = K*w3^2 ./ sqrt((w.^2)*(2*z3*w3)^2 + (w3^2-w.^2).^2);
plot(wr, A3,'-.k','linewidth',2),hold on

w = w4 * wr;
A4 = K*w4^2 ./ sqrt((w.^2)*(2*z4*w4)^2 + (w4^2-w.^2).^2);
plot(wr, A4,'-k','linewidth',2),w = w5 * wr;

w = w5 * wr;
A5 = K*w5^2 ./ sqrt((w.^2)*(2*z5*w5)^2 + (w5^2-w.^2).^2);
plot(wr, A5,':k','linewidth',2),grid, xlabel('w/w_n')

title('Amplitudes')
legend(sprintf('Heave: damping z_3 = %3.2f, nat. frequency w_3 = %3.2f',z3, w3),...
       sprintf('Roll:  damping z_4 = %3.2f, nat. frequency w_4 = %3.2f',z4, w4),...
       sprintf('Pitch: damping z_5 = %3.2f, nat. frequency w_5 = %3.2f',z5, w5))
hold off