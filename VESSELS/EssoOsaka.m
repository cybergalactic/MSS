% Esso Ossaka Nonlinear Maneuvering Mathematical Model
%
% Reference : Abkowitz, Martin (1980). Measurement of Hydrodynamic Characteristics from Ship Maneuvering Trials by System Identification,
%             Transaction of SNAME, 88:283-318
%
% Author:    Lucia Moreira
% Date:      9th February 2004
% Revisions: 20th February 2004
% 
% Esso Osaka (280000 dwt tanker) particulars and normalization variables - Prime-system I
%Loa = 343.00;     % Length Overall (m)
L = 325.00;        % Length Between Perpendiculars (m) 
B = 53.00;         % Breadth Molded (m)
I = 21.73;         % Draft Molded (m) 
Cb = 0.831;        % Block Coefficient 
nabla = 319400;    % Displacement at Trials (ton)
Ar = 119.817;      % Rudder Area (m^2)
xgorig = 10.30;    % Longitudinal CG at Trials: Forward of Midship (m)
rho = 1025;        % Density of Water
g  = 9.81;	       % Acceleration of Gravity

w = 0.2;           % Wake Fraction - typical value 
d = 9.10;          % Propeller Diameter (m)
Ap = pi*(d/2)^2;   % Propeller Area (m^2)
np = 0.68 ;        % Propeller rps
Umax = 16*0.5144;  % Maximum Speed (m/s)

T=1.0e0;
n=2400;
t=0:T:(n-1)*T;

% Initial conditions
u(1) = 8*0.5144;
v(1) = 0;
r(1)= 0;
psi(1) = 0;
uccor(1)= 8*0.5144;
vccor(1)= 0;
x(1) = 0;
y(1) = 0;

% Parameters and hydrodynamic derivatives
coef1 = 0.0352;  % (m - Yvdot)' 
coef2 = 0.00222; % (Iz - Nrdot)' 
coef3 = 0.0266;  % (Xvr + m)'
coef4 = 0.0188;  % (m - Xudot)'
Yv = -0.0261;
Yr = 0.00365;
Nv = -0.0105;
Nr = -0.00480;
Nd = -0.00283;
N0 = -0.00028;
Nvrr = 0.00611;
Nrrv = 0.00611;
Xee = -0.00224;
Xrrvv = -0.00715;
Xvvrr = -0.00715;
Neee = 0.00116;
Yvrr = -0.0450;
Yrrv = -0.0450; 

rend1 = -0.962e-5;
rend2 = -0.446e-5;
rend3 = 0.0309e-5;

m = 0.0181;                % Nondimensional Mass 
morig = (rho/2)*m*(L^2)*I; % Dimensional Mass 
xg = xgorig/L;

Xvv = -0.006;
Xrr = 0.00515;
Y0 = 1.90e-6;
Yd = 0.00508;
Yeee = -0.00185;

Nvdot = 0.00;
Yrdot = 0.00;
%Xudot = 0.05/((rho/2)*(L^3));

Iz = 0.0625*m*I/L;    % Nondimensional (estimated)

ucor = 0.0;           % Current Magnitude
alfa = 0.0;           % Current Direction or Heading

CR = 0.00226;         % Resistance Coefficient of Ship, including Windmilling Propeller - based onvalue of (u/(n.d)) ?????????????????????
S = L*(1.7*I + B*Cb); % Wetted Surface Area based on Denny-Mumford formula
KT = 0.59;            % Estimated - ??????????????????
k = 0.95;             % Estimated - ??????????????????

for i=1:n-1, 

    if i < 20
        delta = 0;
    else
        delta = 35*pi/180;
    end

Thrust = rend1.*u.^2 + rend2.*np.*u + rend3.*np.^2;

uAinf = -(1 - w).*u + sqrt((1 - w)*(1 - w).*u.^2 + (8/pi)*KT*(np*d)^2);                  % Induced Axial Velocity Behind the Propeller Disk

c = sqrt(((Ap/Ar)*((1 - w).*u + k*uAinf)).^2 + ((Ar - Ap)/Ar)*((1 - w)*(1 - w)).*u.^2);  % Weighted Average Flow Speed Over Rudder

c0 = c(1);                                                                               % Value of the Inflow Velocity to Rudder when the Propeller Rotational Speed and 
                                                                                         % Ship Forward Speed are in Equilibrium in Straight-Ahead Motion (u = u0 and n = n0)

e = delta.*(v./c) + ((r.*L)./(2.*c));                                                    % Effective Rudder Angle

Ucap = sqrt(u.^2 + v.^2);
Ucapcor = sqrt(uccor.^2 + vccor.^2);

f1 = rend1*(rho/2)*(L^2).*(u.^2) + rend2*(rho/2)*(L^3)*np.*u + rend3*(rho/2)*(L^4)*(np^2) - CR*(rho/2)*S.*u.^2 ...
   + Xvv*(rho/2)*(L^2).*(v.^2) + Xee*(rho/2)*(L^2).*(c.^2).*(e.^2) + (Xrr + m*xg)*(rho/2)*(L^4).*(r.^2) + ...
   + (coef3)*(rho/2)*(L^3).*v.*r + Xvvrr*(rho/2)*(L^4)*(1./(Ucapcor.^2)).*(v.^2).*(r.^2);

f2 = Y0*(rho/2)*(L^2).*((uAinf./2).^2) + (Yv*(rho/2)*(L^2).*Ucap.*v + Yd.*(c - c0)*(rho/2)*(L^2).*v) ...
   + ((Yr-m.*(u./Umax))*(rho/2)*(L^3).*Ucap.*r - (Yd/2).*(c - c0).*(rho/2)*(L^3).*r) + Yd*(rho/2)*(L^2).*(c.^2).*delta ...       
   + Yrrv*(rho/2)*(L^4).*(1./Ucap).*(r.^2).*v + Yeee*(rho/2)*(L^2).*(c.^2).*(e.^3);

f3 = N0*(rho/2)*(L^3).*((uAinf./2).^2) + (Nv*(rho/2)*(L^3).*Ucap.*v - Nd.*(c - c0)*(rho/2)*(L^3).*v) ...
   + ((Nr-m*xg.*(u./Umax))*(rho/2)*(L^4).*Ucap.*r + (Nd/2).*(c - c0)*(rho/2)*(L^4).*r) + Nd*(rho/2)*(L^3).*(c.^2).*delta ...       
   + Nrrv*(rho/2)*(L^5).*(1./Ucap).*(r.^2).*v + Neee*(rho/2)*(L^3).*(c.^2).*(e.^3);

f4 = (coef1)*(rho/2)*(L^3)*(coef2)*(rho/2)*(L^5) ...
   - (m*xg - Nvdot)*(rho/2)*(L^4)*(m*xg - Yrdot)*(rho/2)*(L^4); 

u(:,i+1) = u(:,i) + T.*(f1(:,i)./((coef4)*(rho/2)*(L^3)));
v(:,i+1) = v(:,i) + T.*(((coef2)*(rho/2)*(L^5).*f2(:,i) - (m*xg - Yrdot)*(rho/2)*(L^4).*f3(:,i))./f4);
r(:,i+1) = r(:,i) + T.*(((coef1)*(rho/2)*(L^3).*f3(:,i) - (m*xg - Nvdot)*(rho/2)*(L^4).*f2(:,i))./f4);
psi(:,i+1) = psi(:,i) + T.*r(:,i);

uccor(:,i+1) = uccor(:,i) + T.*(u(:,i+1) - ucor.*r(:,i)*sin(psi(:,i) - alfa)); 
vccor(:,i+1) = vccor(:,i) + T.*(v(:,i+1) - ucor.*r(:,i)*cos(psi(:,i) - alfa)); 

x(:,i+1) = x(:,i) + T.*(u(:,i).*cos(psi(:,i)) - v(:,i).*sin(psi(:,i)));
y(:,i+1) = y(:,i) + T.*(u(:,i).*sin(psi(:,i)) + v(:,i).*cos(psi(:,i)));

end

udim = u./0.5144;
vdim = v./0.5144;
rdim = r.*180/pi;

figure(1)
subplot(3,1,1)
plot(t,udim)
grid
title('Surge')
ylabel('knots')
subplot(3,1,2)
plot(t,vdim)
grid
title('Sway')
ylabel('knots')
subplot(3,1,3)
plot(t,rdim)
grid
title('Yaw')
xlabel('Time (s)')
ylabel('deg/s')

figure(2)
plot(x,y)
grid
xlabel('X (m)')
ylabel('Y (m)')
title('Simulation of Turning Circle (Max.Rudder, Half Speed)')
