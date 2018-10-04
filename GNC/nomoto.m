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

T = T1+T2-T3;
d1 = [T 1 0];           n1 = K; 
d2 = [T1*T2 T1+T2 1 0]; n2 = K*[T3 1];
[mag1,phase1,w] = bode(n1,d1);
[mag2,phase2]   = bode(n2,d2,w);

if K < 0,   % shift ship phase with 360 deg for course unstable ship
phase1 = phase1-360;
phase2 = phase2-360;
end

clf,subplot(211),semilogx(w,20*log10(mag1),'b'),grid
xlabel('Frequency [rad/s]'),title('Gain [dB]')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(212),semilogx(w,phase1,'b'),grid
xlabel('Frequency [rad/s]'),title('Phase [deg]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')