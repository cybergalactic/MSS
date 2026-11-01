function nomoto(T1,T2,T3,K)
% NOMOTO(T1,T2,T3,K) generates the Bode plots for
%
%              K                       K (1+T3s)
%  H1(s) = ---------       H2(s) = -------------------
%           (1+Ts)s                  s(1+T1s)(1+T2s)
%
% Author:   Thor I. Fossen
% Date:     19th June 2001
% Revisions: 
%   2026-11-01 - Direct frequency response for base MATLAB and GNU Octave.

T = T1 + T2 - T3;
w = logspace(-4,1,500);
s = 1i * w;
H1 = K ./ (s .* (1 + T*s));
H2 = K * (1 + T3*s) ./ (s .* (1 + T1*s) .* (1 + T2*s));
mag1 = abs(H1);
mag2 = abs(H2);
phase1 = unwrap(angle(H1)) * 180/pi;
phase2 = unwrap(angle(H2)) * 180/pi;

if K < 0   % Shift ship phase with 360 deg for course unstable ship
phase1 = phase1 - 360;
phase2 = phase2 - 360;
end

clf,subplot(211),semilogx(w,20*log10(mag1),'b--'),grid
xlabel('Frequency [rad/s]'),title('Gain [dB]')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(212),semilogx(w,phase1,'b--'),grid
xlabel('Frequency [rad/s]'),title('Phase [deg]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')