% exNomoto requires the Matlab control systems toolbox
% The script generates the Nomoto Bode plots of two MSS ships.
%      
% Author:   Thor I. Fossen
% Date:     2001-06-19
% Revisions: 

% Mariner-class cargo ship
T1 = 118; T2 = 7.8; T3 = 18.5; K = 0.185;  

T = T1+T2-T3;
d1 = [T 1 0];           n1 = K; 
d2 = [T1*T2 T1+T2 1 0]; n2 = K*[T3 1];
w =logspace(-4,1,50);
wi = [0.0001 10];
[mag1,phase1] = bode(n1,d1,w);
[mag2,phase2] = bode(n2,d2,w);

if K < 0    % Shift ship phase with 360 deg for course unstable ship
    phase1 = phase1 - 360;
    phase2 = phase2 - 360;
end

clf
subplot(411)
semilogx(w,20*log10(mag1),'b')
grid
ylabel('Gain [dB]'),title('Mariner-class cargo vessel')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(412),semilogx(w,phase1,'b'),grid
ylabel('Phase [deg]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')

% Oil tanker
T1 = -124.1; T2 = 16.4; T3 = 46.0; K = -0.019; 

T = T1+T2-T3;
d1 = [T 1 0];           n1 = K; 
d2 = [T1*T2 T1+T2 1 0]; n2 = K*[T3 1];
[mag1,phase1] = bode(n1,d1,w);
[mag2,phase2] = bode(n2,d2,w);

if K < 0   % Shift ship phase with 360 deg for course unstable ship
    phase1 = phase1 - 360;
    phase2 = phase2 - 360;
end

subplot(413),semilogx(w,20*log10(mag1),'b'),grid
ylabel('Gain [dB]'),title('Tanker')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(414),semilogx(w,phase1,'b'),grid
ylabel('Phase [deg]');xlabel('Frequency [rad/s]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
