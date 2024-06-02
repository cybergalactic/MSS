function [xdot, U, B_delta] = npsauv(x, ui)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [xdot, U, B_delta] = npsauv(x, ui) returns the time derivative of the 
% state vector: x = [u v w p q r xpos ypos zpos phi theta psi
% delta_r delta_s delta_bp delta_bs n]', speed U in m/s (optionally) and, 
% input matrix B_delta (optionally) for an Autonomous Underwater Vehicle 
% (AUV) at the Naval Postgraduate School, Monterrey, USA. The length of the 
% AUV is 5.3 m and the mass is 5443 kg, while the state vector is defined by:
%
%   u:        Surge velocity              (m/s)
%   v:        Sway velocity               (m/s)
%   w:        Heave velocity              (m/s)
%   p:        Roll rate                   (rad/s)
%   q:        Pitch rate                  (rad/s)
%   r:        Yaw rate                    (rad/s)
%   xpos:     Position in x-direction     (m)
%   ypos:     Position in y-direction     (m)
%   zpos:     Position in z-direction     (m)
%   phi:      Roll angle                  (rad)
%   theta:    Pitch angle                 (rad)
%   psi:      Yaw angle                   (rad)
%   delta_r:  Rudder angle                (rad)
%   delta_s:  Stern plane                 (rad)
%   delta_bp: Port bow plane              (rad)
%   delta_bs: Starboard bow plane         (rad)
%   n:        Propeller shaft speed       (rpm)
%
% Input commands:
%
%   ui = [ delta_r_com delta_s_com delta_bp_com delta_bs_com n_com ]'  
%
%   delta_r_com:    Rudder angle command          (rad)
%   delta_s_com:    Stern plane command           (rad)
%   delta_bp_com:   Port bow plane command        (rad)
%   delta_bs_com:   Starboard bow plane command   (rad)
%   n_com:          Propeller shaft speed command (rpm)  
%
% M-file Simulators:
%   SIMnpsauv.m : Script demonstrating 3-D ALOS path-following control.
%
% Simulink Simulators:
%   demoNPSAUV.slx : Simulink model demonstrating PID heading control.
%
% Reference: 
%   A. J. Healey and Lienard, D. (1993). Multivariable Sliding Mode Control 
%     for Autonomous Diving and Steering of Unmanned Underwater Vehicles,
%     IEEE Journal of Ocean Engineering 18(3):327-339.
%
% Author:    Thor I. Fossen
% Date:      2024-06-03
% Revisions: 

% Check of input and state dimensions
if (length(x) ~= 17), error('x-vector must have dimension 17!'); end
if (length(ui) ~= 5), error('u-vector must have dimension 5!'); end

% Dimensional states
u   = x(1);  v     = x(2);  w   = x(3);
p   = x(4);  q     = x(5);  r   = x(6);
phi = x(10); theta = x(11); psi = x(12);
u_actual = x(13:17);

U = sqrt(u^2 + v^2 + w^2);  % speed

% Actuator dynamics and saturation limits
max_ui = zeros(5,1);
max_ui(1) = deg2rad(20);   % Max delta_r   (rad)
max_ui(2) = deg2rad(20);   % Max delta_s   (rad)
max_ui(3) = deg2rad(20);   % Max delta_bp  (rad)
max_ui(4) = deg2rad(20);   % Max delta_bs  (rad)
max_ui(5) = 1500;          % Max propeller speed (rpm)

u_actual = sat(u_actual, max_ui);   % Saturation of control inputs 

T_actuator = 1.0;          % Unified actuator time constant (seconds)
udot = (ui - u_actual) / T_actuator;  % Actuator dynamics

delta_r  = u_actual(1);    % Actual control inputs
delta_s  = u_actual(2);
delta_bp = u_actual(3);
delta_bs = u_actual(4);
n        = u_actual(5) / 60 * 2*pi;

% AUV parameters
L   = 5.3;    g   = 9.81;
xG  = 0;      yG  = 0;      zG = 0.061;
xB  = 0;      yB  = 0;      zB = 0;
W   = 53400;  B   = 53400;
rho = 1025;   mass = W/g;   
Ix  = 2038;   Iy  = 13587;  Iz  = 13587;
Ixy = -13.58; Iyz = -13.58; Ixz = -13.58;
Cdy = 0.5;    Cdz = 0.6;
r_bb = [xB yB zB]';
r_bg = [xG yG zG]';

% Prime-scaling variables (Fossen 2021, Appendix D.2)
Tinv = diag( [1 1 1 L L L] );
r2 = (1/2) * rho * L^2;
r3 = (1/2) * rho * L^3;
r4 = (1/2) * rho * L^4;
r5 = (1/2) * rho * L^5;

% Nondimensional mass and hydrodynamic derivaties
m = mass / r3;

Xpp   =  7.0e-3; Xqq    = -1.5e-2; Xrr   =  4.0e-3; Xpr   =  7.5e-4;
Xudot = -7.6e-3; Xwq    = -2.0e-1; Xvp   = -3.0e-3; Xvr   =  2.0e-2;
Xqds  =  2.5e-2; Xqdb2  = -1.3e-3; Xrdr  = -1.0e-3; Xvv   =  5.3e-2;
Xww   =  1.7e-1; Xvdr   =  1.7e-3; Xwds  =  4.6e-2; Xwdb2 =  0.5e-2;
Xdsds = -1.0e-2; Xdrdr  = -1.0e-2; Xqdsn =  2.0e-3;
Xwdsn =  3.5e-3; Xdsdsn = -1.6e-3;

Ypdot =  1.2e-4; Yrdot  =  1.2e-3; Ypq   =  4.0e-3; Yqr   = -6.5e-3;
Yvdot = -5.5e-2; Yp     =  3.0e-3; Yr    =  3.0e-2; Yvq   =  2.4e-2;
Ywp   =  2.3e-1; Ywr    = -1.9e-2; Yv    = -1.0e-1; Yvw   =  6.8e-2;
Ydr   =  2.7e-2;

Zqdot = -6.8e-3; Zpp    =  1.3e-4; Zpr   =  6.7e-3; Zrr   = -7.4e-3;
Zwdot = -2.4e-1; Zq     = -1.4e-1; Zvp   = -4.8e-2; Zvr   =  4.5e-2;
Zw    = -3.0e-1; Zvv    = -6.8e-2; Zds   = -7.3e-2; Zdb2  = -1.3e-2;
Zqn   = -2.9e-3; Zwn    = -5.1e-3; Zdsn  = -1.0e-2;

Kpdot = -1.0e-3; Krdot  = -3.4e-5; Kpq   = -6.9e-5; Kqr   =  1.7e-2;
Kvdot =  1.2e-4; Kp     = -1.1e-2; Kr    = -8.4e-4; Kvq   = -5.1e-3;
Kwp   = -1.3e-4; Kwr    =  1.4e-2; Kv    =  3.1e-3; Kvw   = -1.9e-1;
Kdb2  =  0;      Kpn    = -5.7e-4; Kprop =  0;

Mqdot = -1.7e-2; Mpp    =  5.3e-5; Mpr   =  5.0e-3; Mrr   =  2.9e-3;
Mwdot = -6.8e-3; Muq    = -6.8e-2; Mvp   =  1.2e-3; Mvr   =  1.7e-2;
Muw   =  1.0e-1; Mvv    = -2.6e-2; Mds   = -4.1e-2; Mdb2  =  3.5e-3;
Mqn   = -1.6e-3; Mwn    = -2.9e-3; Mdsn  = -5.2e-3;

Npdot = -3.4e-5; Nrdot  = -3.4e-3; Npq   = -2.1e-2; Nqr   =  2.7e-3;
Nvdot =  1.2e-3; Np     = -8.4e-4; Nr    = -1.6e-2; Nvq   = -1.0e-2;
Nwp   = -1.7e-2; Nwr    =  7.4e-3; Nv    = -7.4e-3; Nvw   = -2.7e-2;
Ndr   = -1.3e-2; Nprop  =  0;

% Rigid-body and added mass matrices
I_g = [Ix -Ixy -Ixz
       -Ixy Iy -Iyz
       -Ixz -Iyz Iz ];

MRB = [  mass * eye(3)      -mass * Smtrx(r_bg)
         mass * Smtrx(r_bg)  I_g                ];

MA = -[ Xudot 0     0     0     0     0
        0     Yvdot 0     Ypdot 0     Yrdot
        0     0     Zwdot 0     Zqdot 0
        0     Kvdot 0     Kpdot 0     Krdot
        0     0     Mwdot 0     Mqdot 0
        0     Nvdot 0     Npdot 0     Nrdot ];
MA = r3 * Tinv * MA * Tinv;   

M = MRB + MA;

% Control forces and moments
Cd0   = 0.00385;
eps_prop = 1e-10;  % Avoid dividing by zero
prop  = 0.012 * n / (u + eps_prop); 
Xprop = Cd0 * ( abs(prop) * prop - 1 );
Ct    = 0.008 * L^2 * abs(prop) * prop / 2;
Ct1   = 0.008 * L^2 / 2;
epsilon = -1 + sign(n)/sign(u) * ( sqrt(Ct+1)-1) / (sqrt(Ct1+1)-1 );

tau1 = r3 * (Xrdr*u*r*delta_r + (Xqds*delta_s + Xqdb2*delta_bp +...
             Xqdb2*delta_bs)*u*q) +...
       r2 * (Xvdr*u*v*delta_r + (Xwds*delta_s + Xwdb2*delta_bs + ...
             Xwdb2*delta_bp)*u*w + (Xdsds*delta_s^2 + ...
             + Xdrdr*delta_r^2)*u^2) +...
       r3 * Xqdsn*u*q*delta_s*epsilon + r2*(Xwdsn*u*w*delta_s +...
            Xdsdsn*u^2*delta_s^2)*epsilon + r2*u^2*Xprop;
tau2 = r2 * Ydr*u^2*delta_r;
tau3 = r2 * u^2*(Zds*delta_s+Zdb2*delta_bs+Zdb2*delta_bp) + ...
       r3 * Zqn*u*q*epsilon+r2*(Zwn*u*w+Zdsn*u^2*delta_s)*epsilon;
tau4 = r4 * Kpn*u*p*epsilon+r3*u^3*Kprop + ...
       r3 * u^2*(Kdb2*delta_bp+Kdb2*delta_bs);
tau5 = r4 * Mqn*u*q*epsilon+r3*(Mwn*w*u+Mdsn*u^2*delta_s)*epsilon+...
       r3 * u^2*(Mds*delta_s+Mdb2*delta_bp+Mdb2*delta_bs);
tau6 = r3 * u^2*Nprop + r3*u^2*Ndr*delta_r;

tau_control = [tau1 tau2 tau3 tau4 tau5 tau6]';

% 4-DOF optional input matrix B_delta for the control surfaces in sway, 
% heave, pitch, and yaw, defined by
%   tau_control = [ tau1
%                   B_delta * u_delta + b0(ui) ]
%   u_delta = [ delta_r delta_s delta_bp delta_bs ]' 
%   b0 = [ tau1
%          0
%          r3 * Zqn * u * q * epsilon + r2 * Zwn * u * w * epsilon
%          tau4
%          r4 * Mqn * u * q * epsilon + r3 * Mwn * u * w * epsilon
%          r3 * u^2 * Nprop ]
B_delta = u^2 * [
    r2 * Ydr, 0, 0, 0
    0, r2 * Zds + r2 * Zdsn * epsilon, r2 * Zdb2, r2 * Zdb2
    0, r3 * Mdsn * epsilon + r3 * Mds, r3 * Mdb2, r3 * Mdb2
    r3 * Ndr, 0, 0, 0 ];

% Cross-flow drag assuming a block-shaped body and strip theory
dxL = L/10;             % Number of sections
epsilon_drag = 1e-6;    % Avoid dividing by zero
Hx = 0.53;              % Average height of submerged body
Bx = 0.53;              % Average width of submerged body

Cy = 0; Cz = 0; Cm = 0; Cn = 0;
for xL = -L/2:dxL:L/2   % Cross-flow drag integrals

    Ucf = sqrt((v + xL * r)^2 + (w - xL * q)^2) + epsilon_drag;
    drag_term = Cdy * Hx * (v + xL * r)^2 + Cdz * Bx * (w - xL * q)^2;
    
    Cy = Cy + dxL * drag_term * (v + xL * r) / Ucf;
    Cz = Cz + dxL * drag_term * (w - xL * q) / Ucf;
    Cm = Cm + dxL * drag_term * (w + xL * q) / Ucf * xL;
    Cn = Cn + dxL * drag_term * (v + xL * r) / Ucf * xL;

end
tau_crossflow = (rho/2) * [ 0 -Cy Cz 0 -Cm -Cn]';

% Restoring forces and moments
tau_hydrostatic = -gvect(W, B, theta, phi, r_bg, r_bb);

% Hydrodynamic forces and moments
X_h = r2 * ( Xvv*v^2 + Xww*w^2 ) + ...
      r3 * ( (m+Xvr)*v*r + (Xwq-m)*w*q + Xvp*v*p ) + ... 
      r4 * ( (m*xG/L+Xqq)*q^2 + (m*xG/L+Xrr)*r^2 - m*yG/L*p*q + ...
             (Xpr-m*zG/L)*p*r + Xpp*p^2 );

Y_h = r2 * ( Yv*u*v + Yvw*v*w ) + ...
      r3 * ( Yp*u*p + Yr*u*r +Yvq*v*q + Ywp*w*p + Ywr*w*r ) + ...
      r4 * ( Ypq*p*q + Yqr*q*r) - ...
      mass * ( u*r -w*p + xG*p*q - yG*(p^2+r^2) + zG*q*r );

Z_h = r2 * ( Zw*w*u + Zvv*v^2 ) +...
      r3 * ( Zq*u*q + Zvp*v*p + Zvr*v*r ) +...
      r4 * ( Zpp*p^2 + Zpr*p*r + Zrr*r^2 ) + ...
      mass * ( v*p - u*q + xG*p*r + yG*q*r -zG*(p^2+q^2) );

K_h = r3 * ( Kv*u*v + Kvw*v*w ) +...
      r4 * ( Kp*u*p + Kr*u*r + Kvq*v*q + Kwp*w*p + Kwr*w*r ) +...
      r5 * ( Kpq*p*q + Kqr*q*r ) +...
      (Iy-Iz) * q*r - Ixy * p*r - (r^2-q^2) * Iyz + Ixz * p*q - ...
      mass * ( yG*(v*p-u*q) - zG*(u*r - w*p) );

M_h = r3 * ( Muw*u*w + Mvv*v^2 ) +...
      r4 * ( Muq*u*q + Mvp*v*p + Mvr*v*r ) +...
      r5 * ( Mpp*p^2 + Mpr*p*r + Mrr*r^2 ) - ...
      (Iz-Ix) * p*r + Ixy * q*r - Iyz * p*q - (p^2-r^2) * Ixz + ...
      mass * ( xG*(v*p-u*q) - zG*(w*q-v*r) );
 
N_h = r3 * ( Nv*u*v + Nvw*v*w )+...
      r4 * ( Np*u*p + Nr*u*r + Nvq*v*q + Nwp*w*p + Nwr*w*r )+...
      r5 * ( Npq*p*q + Nqr*q*r ) +...
      (Ix-Iy) * p*q + (p^2-q^2) * Ixy + Iyz * p*r - Ixz * q*r -...
      mass * ( xG*(u*r-w*p) - yG*(w*q-v*r) );

tau_hydrodynamic = [ X_h Y_h Z_h K_h M_h N_h ]';

% State derivatives 
xdot = [ M \ ( tau_control + ...
            tau_hydrodynamic + tau_crossflow + tau_hydrostatic );
         eulerang(phi,theta,psi) * x(1:6) 
         udot ];

end