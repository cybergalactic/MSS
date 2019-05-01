% ExResonance computes the responses in heave, roll and pitch for a marine
% craft exposed to a regular wave. The closed-form solution of a linear
% mass-damper-spring system with sinusoidal forcing F sin(w_e t) is used. 
% Harmonic oscillator:
%   ..     .
% m x  + d x + k x = F sin(w_e t),     w_n = sqrt(k/m)
%
% The results are plotted as change in amplitudes (1/Zm) with respect to the 
% frequency ratio w_e/w_n where Zm denotes the impedance.
%
% Author:    Thor I. Fossen
% Date:      2019-05-01
% Revisions: 

load supply;		% load supply data (ShipX)

% compute relative damping rations and natural frequencies
[T,zeta] = DPperiods(vessel);
z3 = zeta(3);   z4 = zeta(4);   z5 = zeta(5);
w3 = 2*pi/T(3); w4 = 2*pi/T(4); w5 = 2*pi/T(5);

wr = 0:0.01:2.5;   % wr = w_e / w_n

figure(gcf)

w = w3 * wr;
Zm3w = sqrt((w.^2)*(2*z3*w3)^2 + (w3^2-w.^2).^2);  % Zm3 * w_e
plot(wr, w3^3./Zm3w,'-.k','linewidth',2),hold on

w = w4 * wr;
Zm4w = sqrt((w.^2)*(2*z4*w4)^2 + (w4^2-w.^2).^2);  % Zm4 * w_e
plot(wr, w4^2./Zm4w,'-k','linewidth',2)

z4_double = 2*z4;
Zm4w = sqrt((w.^2)*(2*z4_double*w4)^2 + (w4^2-w.^2).^2);
plot(wr, w4^2./Zm4w,'-k','linewidth',1)

w = w5 * wr;
Zm5w = sqrt((w.^2)*(2*z5*w5)^2 + (w5^2-w.^2).^2);  % Zm5 * w_e
plot(wr, w5^2./Zm5w,':k','linewidth',2)

title('Amplitude = w_i^2 / (Z_m w_e)'),xlabel('w_e/w_n'), grid
legend(sprintf('Heave: damping z_3 = %3.2f, nat. frequency w_3 = %3.2f',z3, w3),...
       sprintf('Roll:  damping z_4 = %3.2f, nat. frequency w_4 = %3.2f',z4, w4),... 
       sprintf('Roll:  damping z_4 = %3.2f, nat. frequency w_4 = %3.2f',z4_double, w4),...        
       sprintf('Pitch: damping z_5 = %3.2f, nat. frequency w_5 = %3.2f',z5, w5))
hold off 