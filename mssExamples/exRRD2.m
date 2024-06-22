% exRRD2: Rudder-Roll Damping (RRD) system for the Son and Nomoto (1982)
% container ship. Calls the function 'RRDcontainer.m' at the end of the 
% file, which includes wave-induced roll disturbances.
%
% References:  
%   K. Son og K. Nomoto (1982). On the Coupled Motion of Steering and 
%   Rolling of a High Speed Container Ship, Naval Architect of 
%   Ocean Engineering 20:73-83. From J.S.N.A., Japan, Vol. 150, 1981.
%
% Author:    Thor I. Fossen
% Date:      2001-10-31
% Revisions: 
%   2021-04-07 - Added function RRDcontainer.m to the end of the file

conversion % Load conversion factors

Ns = 20000; % No. of samples
h  = 0.05;  % Sample time

% Normalization variables
rho = 1025;                 % Water density (kg/m^3)
L = 175;                    % Length of ship (m)
U0 = 7.3;                   % Service speed (m/s)

%% SHIP MODEL
% Lineari model using nondimensional matrices and states with dimension; 
% see Eq. (D.13) in Fossen (2021, Appendix D): 
%   T * M' * inv(T) * v_dot + (U/L) * T * N' * inv(T) * v + 
%      (U/L)^2 * T * G' * inv(T) * eta = (U^2/L) * T * b' * delta
%   where eta = [x y psi]' and nu = [v p r]'
T    = diag([ 1 1/L 1/L]);
Tinv = diag([ 1 L L ]);

M = [ 0.01497    0.0003525     -0.0002205       
     -0.0002205  0.0000210      0              
      0.0003525  0              0.000875  ];

N = [ 0.012035   0.00522    0  
     -0.000314   0.0000075  0.0000692      
      0.0038436  -0.000213  0.00243    ];

G = [ 0 0.0000704  0
      0 0.0004966  0 
      0 0.0001468  0];

b = [-0.002578 0.0000855 0.00126  ]';

% Ship state-space model
Minv = inv(T*M*Tinv);
A11 = - Minv * (U0/L)*T*N*Tinv;
A12 = - Minv * (U0/L)^2*T*G*Tinv;
B1  =   Minv * (U0^2/L)*T*b;

A = [ A11        A12(:,2:3)
      0   1   0   0   0
      0   0   1   0   0    ]
  
B = [ B1 ; 0 ; 0 ]

[PHI,DELTA] = c2d(A,B,h);

C = [0 1 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     0 0 0 0 1 ];

% LQ controller gains
Q = diag([10000 1000 10 1 ]);
R = 5;
[G1,G2] = lqtracker(A,B,C,Q,R)

% Eigenvalues, damping ratios, and natural frequencies
damp(A)
damp(A+B*G1)

%% SIMULATION OF AUTOPILOT AND RRD SYSTEMS
n_ref = 70; % Desired RPM

% x = [u v r x y psi p phi delta n xw1 xw2]'
x = [U0 0 0 0 0 0 0 0 0 n_ref 0 0]';  
xout = zeros(Ns+1,length(x));  

wf1 = 0.1;          % Low-pass filter frequnecy
wf2 = 0.05;         % High-pass filter frequency

xyawf  = [0 0 0]';  % Filter states
xrollf = [0 0]';
xw     = [0 0]';    % Wave states

%% MAIN LOOP 
for i=1:Ns+1
    
    % State vectors for yaw and roll subsystems 
    xyaw  = [x(2) x(3) x(6)]'; 
    xroll = [x(7) x(8)]';
   
    % Discrete time low-pass filter
    phif1 = exp(-h*wf1);
    xyawf = phif1*xyawf + (1-phif1)*xyaw;

    % Discrete time high-pass filter 
    phif2   = exp(-h*wf2);
    yrollf = xrollf + xroll;
    xrollf = phif2*xrollf + (phif2-1)*xroll;
    
    % Filtered states
    x_fb = [xyawf(1) yrollf(1) xyawf(2) yrollf(2) xyawf(3) ]';  

    psi_d = deg2rad(10);  % Desired heading
    
    % LQ optimal control law:  u_ref = G1 * x + G2 * C * x_d
    if i<6000
         % Course-changing, no roll damping
         u_ref = [G1(1) 0 G1(3) 0 G1(5)]*x_fb + G2*[0 0 0 psi_d]';  
    elseif i>6000 && i<15000
        % Course-keeping, roll damping
        u_ref = G1*x_fb + G2*[0 0 0 psi_d]';                        
    else
        % Course-keeping, no roll damping    
        u_ref = [G1(1) 0 G1(3) 0 G1(5)]*x_fb + G2*[0 0 0 psi_d]';   
    end 
   
    % Save simulation results
    xout(i,:) = x';    

    % Nonlinear ship model with wave disturbance
    [xdot,U] = RRDcontainer(x,[u_ref; n_ref]);
    x = x + h*xdot;  
  
end

%% PLOTS
t = h*(0:1:Ns)';
figure(gcf)
subplot(311)
plot(t,R2D*xout(:,6)); 
hold on
plot(t,R2D*psi_d*ones(Ns+1,1),'linewidth',2); 
hold off; 
title('yaw angle \psi (deg)')
subplot(312),plot(t,R2D*xout(:,8)); title('roll angle \phi (deg)')
subplot(313),plot(t,R2D*xout(:,9)); title('rudder angle \delta (deg)')

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

%% Function [xdot,U] = RRDcontainer(x,ui)
function [xdot,U] = RRDcontainer(x,ui)
% [xdot,U] = RRDcontainer(x,ui) returns the speed U in m/s (optionally) and 
% the time derivative of the state vector: 
%   x = [ u v r x y psi p phi delta n xw1 xw2]'  
% for a container ship L = 175 m.
%
% Similar as container.m except that wave disturbances in roll are added
% and that Ddelta_max = 20 instead of 5, where
%
%  u     : surge velocity          (m/s)
%  v     : sway velocity           (m/s)
%  r     : yaw velocity            (rad/s)
%  x     : position in x-direction (m)
%  y     : position in y-direction (m)
%  psi   : yaw angle               (rad)
%  p     : roll velocity           (rad/s)
%  phi   : roll angle              (rad)
%  delta : actual rudder angle     (rad)
%  n     : actual shaft velocity   (rpm)
%  xw1, xw2 : wave states
%
% The input vector is:
%  ui      = [ delta_c n_c ]'  where
%
%  delta_c = commanded rudder angle   (rad)
%  n_c     = commanded shaft velocity (rpm)  
%
% Reference:  
%   K. Son og K. Nomoto (1982). On the Coupled Motion of Steering and 
%   Rolling of a High Speed Container Ship, Naval Architect of 
%   Ocean Engineering 20:73-83. From J.S.N.A., Japan, Vol. 150, 1981. 

%% LINEAR WAVE MODEL 
Kw     = 0.0001;
w_wave = 1.0;                    
Aw = [ 0               1
      -w_wave^2    -2*0.1*w_wave   ];
Bw = [0 Kw ]';

%[PHIW,DELTAW] = c2d(Aw,Bw,h);
% Check of input and state dimensions
if (length(x) ~= 12),error('x-vector must have dimension 12 !');end
if (length(ui) ~= 2),error('u-vector must have dimension  2 !');end

% Normalization variables
L = 175;                     % Length of ship (m)
U = sqrt(x(1)^2 + x(2)^2);   % Service speed (m/s)

% Check service speed
if U <= 0,error('The ship must have speed greater than zero');end
if x(10) <= 0,error('The propeller rpm must be greater than zero');end

delta_max  = 20;             % Max rudder angle (deg)
Ddelta_max = 20;             % Max rudder rate (deg/s)
n_max      = 160;            % Max shaft velocity (rpm)

% Non-dimensional states and inputs
delta_c = ui(1); 
n_c     = ui(2)/60*L/U;  

u     = x(1)/U;   v   = x(2)/U;  
p     = x(7)*L/U; r   = x(3)*L/U; 
phi   = x(8);     psi = x(6); 
delta = x(9);     n   = x(10)/60*L/U;
 
% Parameters, hydrodynamic derivatives and main dimensions
m  = 0.00792;    mx     = 0.000238;   my = 0.007049;
Ix = 0.0000176;  alphay = 0.05;       lx = 0.0313;
ly = 0.0313;     Ix     = 0.0000176;  Iz = 0.000456;
Jx = 0.0000034;  Jz     = 0.000419;   xG = 0;

B     = 25.40;   dF = 8.00;    g     = 9.81;
dA    = 9.00;    d  = 8.50;    nabla = 21222; 
KM    = 10.39;   KB = 4.6154;  AR    = 33.0376;
Delta = 1.8219;  D  = 6.533;   GM    = 0.3/L;
rho   = 1025;    t  = 0.175;   T     = 0.0005; 
 
W     = rho*g*nabla/(rho*L^2*U^2/2);

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

% Rudder saturation and dynamics
if abs(delta_c) >= delta_max*pi/180,
   delta_c = sign(delta_c)*delta_max*pi/180;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max*pi/180,
   delta_dot = sign(delta_dot)*Ddelta_max*pi/180;
end

% Shaft velocity saturation and dynamics
n_c = n_c*U/L;
n   = n*U/L;
if abs(n_c) >= n_max/60,
   n_c = sign(n_c)*n_max/60;
end

if n > 0.3,Tm=5.65/n;else,Tm=18.83;end        
n_dot = 1/Tm*(n_c-n)*60;

% Calculation of state derivatives
vR     = ga*v + cRr*r + cRrrr*r^3 + cRrrv*r^2*v;
uP     = cos(v)*((1 - wp) + tau*((v + xp*r)^2 + cpv*v + cpr*r));
J     = uP*U/(n*D);
KT     = 0.527 - 0.455*J;
uR     = uP*epsilon*sqrt(1 + 8*kk*KT/(pi*J^2));
alphaR = delta + atan(vR/uR);
FN     = - ((6.13*Delta)/(Delta + 2.25))*(AR/L^2)*(uR^2 + vR^2)*sin(alphaR);
T      = 2*rho*D^4/(U^2*L^2*rho)*KT*n*abs(n);

% Forces and moments
X  = Xuu*u^2 + (1-t)*T + Xvr*v*r + Xvv*v^2 + Xrr*r^2 + Xphiphi*phi^2 + ...
    cRX*FN*sin(delta) + (m + my)*v*r;

Y  = Yv*v + Yr*r + Yp*p + Yphi*phi + Yvvv*v^3 + Yrrr*r^3 + Yvvr*v^2*r + ...
    Yvrr*v*r^2 + Yvvphi*v^2*phi + Yvphiphi*v*phi^2 + Yrrphi*r^2*phi + ...
    Yrphiphi*r*phi^2 + (1 + aH)*FN*cos(delta) - (m + mx)*u*r;

K  = Kv*v + Kr*r + Kp*p + Kphi*phi + Kvvv*v^3 + Krrr*r^3 + Kvvr*v^2*r + ...
    Kvrr*v*r^2 + Kvvphi*v^2*phi + Kvphiphi*v*phi^2 + Krrphi*r^2*phi + ...
    Krphiphi*r*phi^2 - (1 + aH)*zR*FN*cos(delta) + mx*lx*u*r - W*GM*phi;

N  = Nv*v + Nr*r + Np*p + Nphi*phi + Nvvv*v^3 + Nrrr*r^3 + Nvvr*v^2*r + ...
    Nvrr*v*r^2 + Nvvphi*v^2*phi + Nvphiphi*v*phi^2 + Nrrphi*r^2*phi + ...
    Nrphiphi*r*phi^2 + (xR + aH*xH)*FN*cos(delta);
     
K = K + x(11);    % Wave-induced roll moments     

% Dimensional state derivatives  xdot = [ u v r x y psi p phi delta n ]'
detM = m22*m33*m44-m32^2*m44-m42^2*m33;

xdot =[                      X*(U^2/L)/m11
          -((-m33*m44*Y+m32*m44*K+m42*m33*N)/detM)*(U^2/L)
           ((-m42*m33*Y+m32*m42*K+N*m22*m33-N*m32^2)/detM)*(U^2/L^2)
                   (cos(psi)*u-sin(psi)*cos(phi)*v)*U
                   (sin(psi)*u+cos(psi)*cos(phi)*v)*U 
                              cos(phi)*r*(U/L)                
           ((-m32*m44*Y+K*m22*m44-K*m42^2+m32*m42*N)/detM)*(U^2/L^2)
                                p*(U/L)
                              delta_dot 
                                n_dot                
                                Aw*x(11:12) + Bw*randn(1) ];
end
