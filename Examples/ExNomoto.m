% ExNomoto  Generates the Nomoto Bode plots of the two ships in Chapter 2
%           See also: nomoto.m
% Author:   Thor I. Fossen
% Date:     19th June 2001
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org

% Mariner class cargo ship
T1 = 118; T2 = 7.8; T3 = 18.5; K = 0.185;  

T = T1+T2-T3;
d1 = [T 1 0];           n1 = K; 
d2 = [T1*T2 T1+T2 1 0]; n2 = K*[T3 1];
w=logspace(-4,1,50);
wi = [0.0001 10];
[mag1,phase1] = bode(n1,d1,w);
[mag2,phase2] = bode(n2,d2,w);

if K < 0,   % shift ship phase with 360 deg for course unstable ship
phase1 = phase1-360;
phase2 = phase2-360;
end

clf,subplot(411),semilogx(w,20*log10(mag1),'b'),grid
ylabel('Gain [dB]'),title('Mariner class vessel')
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

if K < 0,   % shift ship phase with 360 deg for course unstable ship
phase1 = phase1-360;
phase2 = phase2-360;
end

subplot(413),semilogx(w,20*log10(mag1),'b'),grid
ylabel('Gain [dB]'),title('Tanker')
hold on,semilogx(w,20*log10(mag2),'r'),hold off
legend('1st-order model','2nd-order model')
subplot(414),semilogx(w,phase1,'b'),grid
ylabel('Phase [deg]');xlabel('Frequency [rad/s]')
hold on,semilogx(w,phase2,'r'),hold off
legend('1st-order model','2nd-order model')
