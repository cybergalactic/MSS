function loadcond(vessel)
% loadcond   loadcond(vessel) plots T_roll and T_pitch as a function of 
%            GM_T and GM_L. The heave period T_heave is plotted as a 
%            function of draft T.
%  Inputs:
%     vessel: MSS vessel structure
%
% Author:    Thor I. Fossen
% Date:      2005-09-26
%            2019-04-25 Extended to inlcude heave, roll and pitch 
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

% Vessel data
m   = vessel.main.m;
T = vessel.main.T;
R44 = vessel.main.k44;
R55 = vessel.main.k55;
g = 9.81;
Ix = vessel.MRB(4,4);
Iy = vessel.MRB(5,5);

% Vessel periods
per = DPperiods(vessel);
T3 = per(3);  T4 = per(4);  T5 = per(5);
w3 = 2*pi/T3; w4 = 2*pi/T4; w5 = 2*pi/T5;

A44 = interp1(vessel.freqs,reshape(vessel.A(4,4,:),1,length(vessel.freqs)),w4);
A55 = interp1(vessel.freqs,reshape(vessel.A(5,5,:),1,length(vessel.freqs)),w5);
M44 = Ix + A44;
M55 = Iy + A55;

kappa_roll  = A44/(m*R44^2);
kappa_pitch = A55/(m*R55^2);

T_draft = 0.01:0.1:1.5*T;
T_heave = 2*pi*sqrt(2*T_draft/g);
GMT = 1:0.1:1.2*vessel.main.GM_T;
T_roll  = 2*pi*sqrt(M44./(m*g*GMT));
GML = 10:1:1.2*vessel.main.GM_L;
T_pitch  = 2*pi*sqrt(M55./(m*g*GML));

% Plots
figure(100)
subplot(311)
plot(T_draft,T_heave,'-','linewidth',2)
hold on
plot(T,T3,'r*','MarkerSize',7,'linewidth',2)
plot([0 max(T_draft)],[T3 T3],'k-')
hold off
legend(sprintf('T_{heave} = sqrt( 2*T/g )'),...
    sprintf('Design: T = %3.2f m, T = %3.2f s',T,T3))
xlabel('T (m)')
ylabel('T_3 (s)')
title('Heave period T_3 (s) as a function of T (m)')

subplot(312)
plot(GMT,T_roll,'-','linewidth',2)
hold on
plot(vessel.main.GM_T,T4,'r*','MarkerSize',7,'linewidth',2)
plot([0 max(GMT)],[T4 T4],'k-')
hold off
legend(sprintf('k_{roll} = %3.2f',kappa_roll),...
    sprintf('Design: GM_T = %3.2f m, T_4 = %3.2f s',vessel.main.GM_T,T4))
xlabel('GM_T (m)')
ylabel('T_4 (s)')
title('Roll period T_4 (s) as a function of GM_T (m)')

subplot(313)
plot(GML,T_pitch,'-','linewidth',2)
hold on
plot(vessel.main.GM_L,T5,'r*','MarkerSize',7,'linewidth',2)
plot([0 max(GML)],[T5 T5],'k-')
hold off
legend(sprintf('k_{pitch} = %3.2f',kappa_pitch),...
    sprintf('Design: GM_L = %3.2f m, T_5 = %3.2f s',vessel.main.GM_L,T5))
xlabel('GM_L (m)')
ylabel('T_5 (s)')
title('Pitch period T_5 (s) as a function of GM_L (m)')
