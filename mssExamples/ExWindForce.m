% ExWindForce  is compatible with MATLAB and GNU Octave (www.octave.org).
% The script plots the wind coefficients using the Blendermann (1994) and
% Isherwood (1972) formulas for merchant ships.
% 
% Dependencies: 
%            blendermann94
%            isherwood72
%
% Author:    Thor I. Fossen
% Date:      2001-06-06 
% Revisions: 2004-06-06 - New values for A_T and A_L.
%            2008-11-20 - Updated to include several new ships.

clear CX CY CK CN

gam_deg = (0:10:180)';
gamma_r = deg2rad(gam_deg);
V_r = 20;

%%  CASE 1 - Blendermann (1994) 
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
plot(gam_deg, CX, 'b-', gam_deg, CX, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CX)-0.2 max(CX+0.2)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CX','Location','Best')

subplot(222)
plot(gam_deg, CY, 'b-', gam_deg, CY, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CY)-0.1 max(CY+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CY','Location','Best')

subplot(223)
plot(gam_deg, CK, 'b-', gam_deg, CK, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CK)-0.1 max(CK+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CK','Location','Best')

subplot(224)
plot(gam_deg, CN, 'b-', gam_deg, CN, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CN)-0.05 max(CN+0.05)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CN','Location','Best')

if ~isoctave; sgtitle('Blendermann (1994) Wind Coefficients'); end

set(findall(gcf,'type','legend'),'FontSize',14);


%%  CASE 2 - Isherwood (1972)
clear CX CY CK CN

Loa	      =  100;  % Length overall (m)
B	      =  30;   % Beam (m)
ALw       = 900;   % Lateral projected area (m^2)
AFw       = 300;   % Transverse projected area (m^2)
A_SS      = 100;   % Lateral projected area of superstructure (m^2)
S	      = 100;   % Length of perimeter of lateral projection of model (m)
		           % Excluding waterline and slender bodies such as masts and ventilators (m)
C	      = 50;    % Distance from bow of centroid of lateral projected area (m)
M	      = 2;     % Number of distinct groups of masts or king posts seen in lateral
		           % Projection; king posts close against the bridge front are not included

[tau_wind,CX,CY,CN] = isherwood72(gamma_r,V_r,Loa,B,ALw,AFw,A_SS,S,C,M);
CK = 0 * CX;

% Plots
XDATA = [0 30 60 90 120 150 180];

figure(2)
subplot(221)
plot(gam_deg, CX, 'b-', gam_deg, CX, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CX)-0.2 max(CX+0.2)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CX','Location','Best')

subplot(222)
plot(gam_deg, CY, 'b-', gam_deg, CY, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CY)-0.1 max(CY+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CY','Location','Best')

subplot(223)
plot(gam_deg, CK, 'b-', gam_deg, CK, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CK)-0.1 max(CK+0.1)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CK','Location','Best')

subplot(224)
plot(gam_deg, CN, 'b-', gam_deg, CN, 'r*', 'LineWidth', 2);
xlabel('Angle of wind \gamma_w (deg) relative bow','FontSize',11); grid on;
axis([0  max(XDATA) min(CN)-0.05 max(CN+0.05)]);
set(gca,'XTick',XDATA);
set(gca,'FontSize',12)
legend('CN','Location','Best')

if ~isoctave; sgtitle('Isherwood (1972) Wind Coefficients'); end; 

set(findall(gcf,'type','legend'),'FontSize',14);