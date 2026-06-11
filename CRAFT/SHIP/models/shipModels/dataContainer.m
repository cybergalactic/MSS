clc;
%MCC_Container
% Parameters of container ship L = 175 m, where

% Reference:  Son og Nomoto (1982). On the Coupled Motion of Steering and 
%             Rolling of a High Speed Container Ship, Naval Architect of Ocean Engineering,
%             20: 73-83. From J.S.N.A. , Japan, Vol. 150, 1981.
% 
% Author:    Trygve Lauvdal
% Date:      12th May 1994
% Revisions: 18th July 2001 (Thor I. Fossen): added output U, changed order of x-vector
%            20th July 2001 (Thor I. Fossen): changed my = 0.000238 to my = 0.007049

% Parameters, hydrodynamic derivatives and main dimensions
m  = 0.00792;    mx     = 0.000238;   my = 0.007049;
Ix = 0.0000176;  alphay = 0.05;       lx = 0.0313;
ly = 0.0313;     Ix     = 0.0000176;  Iz = 0.000456;
Jx = 0.0000034;  Jz     = 0.000419;   xG = 0;

L = 175;         B  = 25.40;   dF = 8.00;        g = 9.81;
dA    = 9.00;    d  = 8.50;    nabla = 21222; 
KM    = 10.39;   KB = 4.6154;  AR    = 33.0376;
Delta = 1.8219;  D  = 6.533;   GM    = 0.3/L;
rho   = 1025;    t  = 0.175;   T     = 0.0005; 
 
Xuu      = -0.0004226;  Xvr    = -0.00311;    Xrr      = 0.00020; 
Xphiphi  = -0.00020;    Xvv    = -0.00386;

Kv       =  0.0003026;  Kr     = -0.000063;   Kp       = -0.0000075; 
Kphi     = -0.000021;   Kvvv   =  0.002843;   Krrr     = -0.0000462; 
Kvvr     = -0.000588;   Kvrr   =  0.0010565;  Kvvphi   = -0.0012012; 
Kvphiphi = -0.0000793;  Krrphi = -0.000243;   Krphiphi =  0.00003569;

Yv       = -0.0116;     Yr     =  0.00242;    Yp       =  0; 
Yphi     = -0.000063;   Yvvv   = -0.109;      Yrrr     =  0.00177; 
Yvvr     =  0.0214;     Yvrr   = -0.0405;     Yvvphi   =  0.04605;
Yvphiphi =  0.00304;    Yrrphi =  0.009325;   Yrphiphi = -0.001368;

Nv       = -0.0038545;  Nr     = -0.00222;    Np       =  0.000213; 
Nphi     = -0.0001424;  Nvvv   =  0.001492;   Nrrr     = -0.00229; 
Nvvr     = -0.0424;     Nvrr   =  0.00156;    Nvvphi   = -0.019058; 
Nvphiphi = -0.0053766;  Nrrphi = -0.0038592;  Nrphiphi =  0.0024195;

kk     =  0.631;  epsilon =  0.921;  xR    = -0.5;
wp     =  0.184;  tau     =  1.09;   xp    = -0.526; 
cpv    =  0.0;    cpr     =  0.0;    ga    =  0.088; 
cRr    = -0.156;  cRrrr   = -0.275;  cRrrv =  1.96; 
cRX    =  0.71;   aH      =  0.237;  zR    =  0.033;
xH     = -0.48;  

% Masses and moments of inertia
m11 = (m+mx);
m22 = (m+my);
m32 = -my*ly;
m42 = my*alphay;
m33 = (Ix+Jx);
m44 = (Iz+Jz);
detM = m22*m33*m44-m32^2*m44-m42^2*m33;