function [xdot, U, M, B_delta] = npsauv(x, ui, Vc, betaVc, w_c)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [xdot, U, M, B_delta] = npsauv(x, ui) returns the time derivative of the 
% state vector: x = [u v w p q r xpos ypos zpos phi theta psi delta_r 
% delta_s delta_bp delta_bs n]', speed U in m/s (optionally), the 6x6 mass 
% matrix M (optionally), and the 2x4 input matrix B_delta (optionally) in 
% pitch and yaw for an Autonomous Underwater Vehicle (AUV) at the Naval 
% Postgraduate School (Healey and Lienard (1993). The length of the AUV is 
% 5.3 m and the mass is 5443 kg, while the state vector is defined by:
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
%   ui = [ delta_r_com delta_s_com delta_bp_com delta_bs_com n_com ]'  
%
%   delta_r_com:    Rudder angle command          (rad)
%   delta_s_com:    Stern plane command           (rad)
%   delta_bp_com:   Port bow plane command        (rad)
%   delta_bs_com:   Starboard bow plane command   (rad)
%   n_com:          Propeller shaft speed command (rpm)  
%
% The arguments Vc (m/s), betaVc (rad), w_c (m/s) are optional arguments for 
% ocean currents
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), w_c ]  
%
% Example usage:
%   [~,~,M,B_delta] = npsauv() : Return the 6x6 mass matrix M and 2x4 
%                                input matrix B_delta
%   [xdot, U] = npsauv(x,ui)   : Return xdot and U, no ocean currents
%   xdot = npsauv(x,ui)        : Return xdot, no ocean currents
%   xdot = npsauv(x,ui,Vc,betaVc,alphaVc,w_c) : Return xdot, 3-D currents
%   xdot = npsauv(x,ui,Vc,betaVc,alphaVc) : Return xdot, horizontal currents
%
% M-file Simulators:
%   SIMnpsauv.m : Script demonstrating MIMO PID heading and deph control
%                 as well as 3-D ALOS path-following control.
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
% Revisions: 2025-06-05 : Corrected the sign of Cz (crossflow drag).

if nargin == 0, x=zeros(17,1); ui=zeros(5,1); Vc=0; betaVc=0; w_c=0; end 
if (nargin == 2), Vc=0; betaVc=0; w_c=0; end   % no ocean currents
if (nargin == 4), w_c=0; end                   % no vertical ocean currents

if (length(x) ~= 17), error('x-vector must have dimension 17!'); end
if (length(ui) ~= 5), error('u-vector must have dimension 5!'); end

% Dimensional states
u   = x(1);  v     = x(2);  w   = x(3);
p   = x(4);  q     = x(5);  r   = x(6);
phi = x(10); theta = x(11); psi = x(12);
u_actual = x(13:17);

% Speed
U = sqrt(u^2 + v^2 + w^2); 

% Ocean currents expressed in BODY, nu_c = [u_c v_c w_c 0 0 0]'
u_c = Vc * cos( betaVc - psi );                               
v_c = Vc * sin( betaVc - psi );   

u_r = u - u_c;                      % Relative velocities
v_r = v - v_c;
w_r = w - w_c;

Dnu_c = [ v_c*r -u_c*r 0 0 0 0 ]';    % Time derivative of nu_c

% Actuator dynamics and saturation limits
max_u = [deg2rad(20); deg2rad(20); deg2rad(20); deg2rad(20); 1500];
u_actual = sat(u_actual, max_u);   % Saturation of control inputs 
T_actuator = 0.1;                  % Unified actuator time constant in seconds
udot = (ui-u_actual) / T_actuator; % First-order actuator dynamics

delta_r  = u_actual(1);            % Actual control inputs
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

% Rigid-body and added mass matrices (Fossen 2021, Chapter 3)
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
epsilon = -1 + sign(n)/sign(u + eps_prop) * ...
    ( sqrt(Ct+1)-1 ) / ( sqrt(Ct1+1)-1 + eps_prop );

tau1 = r3 * (Xrdr*u_r*r*delta_r + (Xqds*delta_s + Xqdb2*delta_bp +...
             Xqdb2*delta_bs)*u_r*q) +...
       r2 * (Xvdr*u_r*v_r*delta_r + (Xwds*delta_s + Xwdb2*delta_bs + ...
             Xwdb2*delta_bp)*u_r*w_r + (Xdsds*delta_s^2 + ...
             + Xdrdr*delta_r^2)*u_r^2) +...
       r3 * Xqdsn*u_r*q*delta_s*epsilon + r2*(Xwdsn*u_r*w_r*delta_s +...
            Xdsdsn*u_r^2*delta_s^2)*epsilon + r2*u_r^2*Xprop;
tau2 = r2 * Ydr*u_r^2*delta_r;
tau3 = r2 * u_r^2*(Zds*delta_s+Zdb2*delta_bs+Zdb2*delta_bp) + ...
       r3 * Zqn*u_r*q*epsilon+r2*(Zwn*u_r*w_r+Zdsn*u_r^2*delta_s)*epsilon;
tau4 = r4 * Kpn*u_r*p*epsilon+r3*u_r^3*Kprop + ...
       r3 * u_r^2*(Kdb2*delta_bp+Kdb2*delta_bs);
tau5 = r4 * Mqn*u_r*q*epsilon+r3*(Mwn*w_r*u_r+Mdsn*u_r^2*delta_s)*epsilon+...
       r3 * u_r^2*(Mds*delta_s+Mdb2*delta_bp+Mdb2*delta_bs);
tau6 = r3 * u_r^2*Nprop + r3*u_r^2*Ndr*delta_r;

tau_control = [tau1; tau2; tau3; tau4; tau5; tau6];

% 2-DOF constant input matrix B_delta for the control surfaces in pitch and 
% yaw, accessable by: [~,~,M, B_delta] = npsauv()
if nargin == 0
    %   [tau5; tau6] = u^2 * B_delta * u_delta + b0
    %   u_delta = [ delta_r delta_s delta_bp delta_bs ]' 
    %   b0 = [ r4 * Mqn * u * q * epsilon + r3 * Mwn * u * w * epsilon
    %          r3 * u^2 * Nprop ]
    B_delta = [ 0, r3 * Mdsn * epsilon + r3 * Mds, r3 * Mdb2, r3 * Mdb2
                r3 * Ndr, 0, 0, 0 ];
else
    B_delta = [];
end

% Cross-flow drag assuming a block-shaped body and strip theory
dxL = L/10;             % Number of sections
epsilon_drag = 1e-6;    % Avoid dividing by zero
Hx = 0.53;              % Average height of submerged body
Bx = 0.53;              % Average width of submerged body

Cy = 0; Cz = 0; Cm = 0; Cn = 0;
for xL = -L/2:dxL:L/2   % Cross-flow drag integrals

    Ucf = sqrt((v_r + xL * r)^2 + (w_r - xL * q)^2) + epsilon_drag;
    drag_term = Cdy * Hx * (v_r + xL * r)^2 + Cdz * Bx * (w_r - xL * q)^2;
    
    Cy = Cy + dxL * drag_term * (v_r + xL * r) / Ucf;
    Cz = Cz + dxL * drag_term * (w_r - xL * q) / Ucf;
    Cm = Cm + dxL * drag_term * (w_r + xL * q) / Ucf * xL;
    Cn = Cn + dxL * drag_term * (v_r + xL * r) / Ucf * xL;

end
tau_crossflow = (rho/2) * [ 0 -Cy -Cz 0 -Cm -Cn]';

% Restoring forces and moments
tau_hydrostatic = -gvect(W, B, theta, phi, r_bg, r_bb);

% Hydrodynamic forces and moments
X_h = r2 * ( Xvv*v_r^2 + Xww*w_r^2 ) + ...
      r3 * ( (m+Xvr)*v_r*r + (Xwq-m)*w_r*q + Xvp*v_r*p ) + ... 
      r4 * ( (m*xG/L+Xqq)*q^2 + (m*xG/L+Xrr)*r^2 - m*yG/L*p*q + ...
             (Xpr-m*zG/L)*p*r + Xpp*p^2 );

Y_h = r2 * ( Yv*u_r*v_r + Yvw*v_r*w_r ) + ...
      r3 * ( Yp*u_r*p + Yr*u_r*r +Yvq*v_r*q + Ywp*w_r*p + Ywr*w_r*r ) + ...
      r4 * ( Ypq*p*q + Yqr*q*r) - ...
      mass * ( u_r*r - w_r*p + xG*p*q - yG*(p^2+r^2) + zG*q*r );

Z_h = r2 * ( Zw*w_r*u_r + Zvv*v_r^2 ) +...
      r3 * ( Zq*u_r*q + Zvp*v_r*p + Zvr*v_r*r ) +...
      r4 * ( Zpp*p^2 + Zpr*p*r + Zrr*r^2 ) + ...
      mass * ( v_r*p - u_r*q + xG*p*r + yG*q*r -zG*(p^2+q^2) );

K_h = r3 * ( Kv*u_r*v_r + Kvw*v_r*w_r ) +...
      r4 * ( Kp*u_r*p + Kr*u_r*r + Kvq*v_r*q + Kwp*w_r*p + Kwr*w_r*r ) +...
      r5 * ( Kpq*p*q + Kqr*q*r ) +...
      (Iy-Iz) * q*r - Ixy * p*r - (r^2-q^2) * Iyz + Ixz * p*q - ...
      mass * ( yG*(v_r*p-u_r*q) - zG*(u_r*r - w_r*p) );

M_h = r3 * ( Muw*u_r*w_r + Mvv*v_r^2 ) +...
      r4 * ( Muq*u_r*q + Mvp*v_r*p + Mvr*v_r*r ) +...
      r5 * ( Mpp*p^2 + Mpr*p*r + Mrr*r^2 ) - ...
      (Iz-Ix) * p*r + Ixy * q*r - Iyz * p*q - (p^2-r^2) * Ixz + ...
      mass * ( xG*(v_r*p-u_r*q) - zG*(w_r*q-v_r*r) );
 
N_h = r3 * ( Nv*u_r*v_r + Nvw*v_r*w_r )+...
      r4 * ( Np*u_r*p + Nr*u_r*r + Nvq*v_r*q + Nwp*w_r*p + Nwr*w_r*r )+...
      r5 * ( Npq*p*q + Nqr*q*r ) +...
      (Ix-Iy) * p*q + (p^2-q^2) * Ixy + Iyz * p*r - Ixz * q*r -...
      mass * ( xG*(u_r*r-w_r*p) - yG*(w_r*q-v_r*r) );

tau_hydrodynamic = [ X_h Y_h Z_h K_h M_h N_h ]';

% Kinematic transformation matrix - Euler angles
J = eulerang(phi,theta,psi);

% State derivatives - generalized velocity, position, and actuator states
xdot = [ Dnu_c + M \ ( tau_control + ...
            tau_hydrodynamic + tau_crossflow + tau_hydrostatic );
         J * x(1:6) 
         udot ];

end