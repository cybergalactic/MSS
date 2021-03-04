% ExWindForce  plotting the wind coefficients using Isherwoods (1972) formulas
%
% Author:    Thor I. Fossen
% Date:      10 September 2001
% Revisions: 06.06.2004: new values for A_T and A_L
%            20.11.2008: updated to include several new cases

clear CX CY CK CN

gam_deg = (0:10:180)';
gamma_r = (pi/180)*gam_deg;
V_r = 20;

% ****************************************************************
%%  CASE 1 Blendermann (1994) 
% ****************************************************************
CDt       = 0.85;
CDl_AF{1} = 0.55;     % gamma_r = 0
CDl_AF{2} = 0.65;     % gamma_r = pi
delta     = 0.60;
kappa     = 1.4;
AFw       = 160.7;
ALw       = 434.8;
Loa       = 55.0;
sL        = 1.48;
sH        = 5.1;
vessel_no = 13;

[tau_w,CX,CY,CK,CN] = blendermann94(gamma_r,V_r,AFw,ALw,sH,sL,Loa,vessel_no);

% Plots
XDATA = [0 30 60 90 120 150 180];

figure(1)
subplot(221)
plot(gam_deg,CX,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CX)-0.2 max(CX+0.2)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CX','Location','Best')
title('Wind Coefficients')

subplot(222)
plot(gam_deg,CY,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CY)-0.1 max(CY+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CY','Location','Best')

subplot(223)
plot(gam_deg,CK,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CK)-0.1 max(CK+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CK','Location','Best')

subplot(224)
plot(gam_deg,CN,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CN)-0.05 max(CN+0.05)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CN','Location','Best')


% ****************************************************************
%%  CASE 2 Isherwood (1972)
% ****************************************************************
clear CX CY CK CN

Loa	      =  100;  % length overall (m)
B	      =  30;   % beam (m)
ALw       = 900;   % lateral projected area (m^2)
AFw       = 300;   % transverse projected area (m^2)
A_SS      = 100;   % lateral projected area of superstructure (m^2)
S	      = 100;   % length of perimeter of lateral projection of model (m)
		           % excluding waterline and slender bodies such as masts and ventilators (m)
C	      = 50;    % distance from bow of centroid of lateral projected area (m)
M	      = 2;     % number of distinct groups of masts or king posts seen in lateral
		           % projection; king posts close against the bridge front are not included

[tau_wind,CX,CY,CN] = isherwood72(gamma_r,V_r,Loa,B,ALw,AFw,A_SS,S,C,M);
CK = 0*CX;

% Plots
XDATA = [0 30 60 90 120 150 180];

figure(2)
subplot(221)
plot(gam_deg,CX,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CX)-0.2 max(CX+0.2)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CX','Location','Best')
title('Wind Coefficients')

subplot(222)
plot(gam_deg,CY,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CY)-0.1 max(CY+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CY','Location','Best')

subplot(223)
plot(gam_deg,CK,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CK)-0.1 max(CK+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CK','Location','Best')

subplot(224)
plot(gam_deg,CN,'b*-','linewidth',2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CN)-0.05 max(CN+0.05)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CN','Location','Best')