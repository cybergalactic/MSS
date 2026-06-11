function [xdot,U] = tanker(x,ui)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [xdot,U] = tanker(x,ui) returns the speed U in m/s (optionally) and the 
% time derivative of the state vector: x = [ u v r x y psi delta n ]'  for
% a the Esso 190,000-dwt tanker L = 304.8 m (Berlekom and Goddhard 1972, 
% Appendix A) where
%
%   u:     surge velocity, must be positive (m/s) - design speed u = 8.23 m/s
%   v:     sway velocity (m/s)
%   r:     yaw velocity (rad/s)
%   x:     position in x-direction (m)
%   y:     position in y-direction (m)
%   psi:   yaw angle (rad)
%   delta: actual rudder angle (rad)
%   n:     actual shaft velocity (rpm)  - nominal propeller speed is 80 rpm
% 
% The input vector is
%
%   ui = [ delta_c  n_c h ]'  where
%
%   delta_c: commanded rudder angle (rad)
%   n_c:     commanded shaft velocity (rpm)
%   h:       water depth, must be larger than draft (m) - draft is 18.46 m
%
% Reference: 
%   W. B. Van Berlekom and T. A. and Goddard (1972). Maneuvering of Large
%     Tankers, Transaction of SNAME, 80:264-298
%
% Author:    Trygve Lauvdal
% Date:      1994-05-12
% Revisions: 
%   2001-07-20 : Added speed output U, changed order of x-vector.
%   2005-05-02 : Changed the incorrect expression 
%                c = sqrt(cun^2*u*n + cnn^2*n^2) to 
%                c = sqrt(cun*u*n + cnn*n^2)
%   2024-04-19 : Added compability to GNU Octave.

if (length(x)  ~= 8),error('x-vector must have dimension 8!'); end
if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end

% Normalization variables
L   =  304.8;          % Length of ship (m)
g   =  9.8;            % Acceleration of gravity (m/s^2)

% Dimensional states and inputs
u     = x(1);    
v     = x(2); 
r     = x(3);
psi   = x(6); 
delta = x(7);
n     = x(8) / 60; % rps
U     = sqrt(x(1)^2 + x(2)^2);

delta_c = ui(1); 
n_c     = ui(2) / 60;  % rps
h       = ui(3);

% Parameters, hydrodynamic derivatives and main dimensions
delta_max  = 10;       % Max rudder angle      (deg)
Ddelta_max = 2.33;     % Max rudder derivative (deg/s)
n_max      = 80;       % Max shaft speed   (rpm)

t   =  0.22;
Tm  =  50;
T   =  18.46;

cun =  0.605;  
cnn =  38.2;

Tuu = -0.00695;
Tun = -0.00063;
Tnn =  0.0000354;

m11 =  1.050;          % 1 - Xudot
m22 =  2.020;          % 1 - Yvdot
m33 =  0.1232;         % kz^2 - Nrdot

d11 =  2.020;          % 1 + Xvr
d22 = -0.752;          % Yur - 1
d33 = -0.231;          % Nur - xG 

Xuuz   = -0.0061;   YT     =  0.04;   NT      = -0.02;
Xuu    = -0.0377;   Yvv    = -2.400;  Nvr     = -0.300;
Xvv    =  0.3;      Yuv    = -1.205;  Nuv     = -0.451;   
Xudotz = -0.05;     Yvdotz = -0.387;  Nrdotz  = -0.0045;
Xuuz   = -0.0061;   Yurz   =  0.182;  Nurz    = -0.047;
Xvrz   =  0.387;    Yvvz   = -1.5;    Nvrz    = -0.120;
Xccdd  = -0.093;    Yuvz   =  0;      Nuvz    = -0.241;
Xccbd  =  0.152;    Yccd   =  0.208;  Nccd    = -0.098;
Xvvzz  =  0.0125;   Yccbbd = -2.16;   Nccbbd  =  0.688;
                    Yccbbdz= -0.191;  Nccbbdz =  0.344;

% Additional terms in shallow water
z = T / (h - T);
if h < 18.5 
    error('The depth must be larger than the draft (18.5 m)'); 
end
if z >= 0.8
    Yuvz = -0.85 * (1 - 0.8/z);
end 

% Rudder saturation and dynamics
if abs(delta_c) >= deg2rad(delta_max)
   delta_c = sign(delta_c) * deg2rad(delta_max);
end
delta_dot = delta_c - delta;
if abs(delta_dot) >= deg2rad(Ddelta_max)
   delta_dot = sign(delta_dot) * deg2rad(Ddelta_max);
end

% Shaft saturation and dynamics
if abs(n_c) >= n_max/60
   n_c = sign(n_c) * n_max/60;
end                 

n_dot = 1/Tm*(n_c-n)*60;

% Forces and moments
if u <= 0
    error('u must be larger than zero'); 
end

beta = atan(v / u);
gT   = (1/L*Tuu*u^2 + Tun*u*n + L*Tnn*abs(n)*n);
c    = sqrt(cun*u*n + cnn*n^2);

gX   = 1/L*(Xuu*u^2 + L*d11*v*r + Xvv*v^2 + Xccdd*abs(c)*c*delta^2 ...
     + Xccbd*abs(c)*c*beta*delta + L*gT*(1-t) + Xuuz*u^2*z ...
     + L*Xvrz*v*r*z + Xvvzz*v^2*z^2);

gY   = 1/L*(Yuv*u*v + Yvv*abs(v)*v + Yccd*abs(c)*c*delta + L*d22*u*r ...
     + Yccbbd*abs(c)*c*abs(beta)*beta*abs(delta) + YT*gT*L ...
     + L*Yurz*u*r*z + Yuvz*u*v*z + Yvvz*abs(v)*v*z  ...
     + Yccbbdz*abs(c)*c*abs(beta)*beta*abs(delta)*z);     

gLN  = Nuv*u*v + L*Nvr*abs(v)*r + Nccd*abs(c)*c*delta +L*d33*u*r ...
     + Nccbbd*abs(c)*c*abs(beta)*beta*abs(delta) + L*NT*gT ...
     + L*Nurz*u*r*z + Nuvz*u*v*z + L*Nvrz*abs(v)*r*z ...
     + Nccbbdz*abs(c)*c*abs(beta)*beta*abs(delta)*z;

m11 = (m11 - Xudotz*z);
m22 = (m22 - Yvdotz*z);
m33 = (m33 - Nrdotz*z);

% Dimensional state derivative
xdot = [  gX/m11
          gY/m22
          gLN/(L^2*m33)
          cos(psi) * u - sin(psi) * v
          sin(psi) * u + cos(psi) * v
          r
          delta_dot
          n_dot    ];