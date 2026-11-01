% exNomoto uses base MATLAB or GNU Octave; no control toolbox is required.
% The script generates the Nomoto Bode plots of two MSS ships.
%      
% Author:   Thor I. Fossen
% Date:     2001-06-19
% Revisions: 
%   2026-11-01 - Direct frequency response with continuous phase branches.

% Mariner-class cargo ship
T1 = 118; T2 = 7.8; T3 = 18.5; K = 0.185;  

T = T1 + T2 - T3;
w = logspace(-4,1,500);
s = 1i * w;
H1 = K ./ (s .* (1 + T*s));
H2 = K * (1 + T3*s) ./ (s .* (1 + T1*s) .* (1 + T2*s));
mag1 = abs(H1);
mag2 = abs(H2);
phase1 = unwrap(angle(H1)) * 180/pi;
phase2 = unwrap(angle(H2)) * 180/pi;

if K < 0    % Shift ship phase with 360 deg for course unstable ship
    phase1 = phase1 - 360;
    phase2 = phase2 - 360;
end

clf
subplot(411)
semilogx(w,20*log10(mag1),'b--')
grid
ylabel('Gain [dB]'),title('Mariner-class cargo vessel')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(412),semilogx(w,phase1,'b--'),grid
ylabel('Phase [deg]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')

% Oil tanker
T1 = -124.1; T2 = 16.4; T3 = 46.0; K = -0.019; 

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

subplot(413),semilogx(w,20*log10(mag1),'b--'),grid
ylabel('Gain [dB]'),title('Tanker')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(414),semilogx(w,phase1,'b--'),grid
ylabel('Phase [deg]');xlabel('Frequency [rad/s]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
