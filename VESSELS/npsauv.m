function [xdot,U] = npsauv(x,ui)
% [xdot,U] = NPSAUV(x,ui) returns the speed U in m/s (optionally) and the 
% time derivative of the state vector: x = [ u v w p q r x y z phi theta psi ]' for
% an Autnomous Underwater Vehicle (AUV) at the Naval Postgraduate School, Monterrey.
% The length of the AUV is L = 5.3 m, while the state vector is defined as:
%
% u     = surge velocity          (m/s)
% v     = sway velocity           (m/s)
% w     = heave velocity          (m/s)
% p     = roll velocity           (rad/s)
% q     = pitch velocity          (rad/s)
% r     = yaw velocity            (rad/s)
% xpos  = position in x-direction (m)
% ypos  = position in y-direction (m)
% zpos  = position in z-direction (m)
% phi   = roll angle              (rad)
% theta = pitch angle             (rad)
% psi   = yaw angle               (rad)
%
% The input vector is :
%
% ui       = [ delta_r delta_s delta_b delta_bp delta_bs n ]'  where
%
% delta_r  = rudder angle                   (rad)
% delta_s  = port and starboard stern plane (rad)
% delta_b  = top and bottom bow plane       (rad)
% delta_bp = port bow plane                 (rad)
% delta_bs = starboard bow plane            (rad)
% n        = propeller shaft speed          (rpm)  
%
% Reference : Healey, A.J. and Lienard, D. (1993). Multivariable Sliding Mode Control 
%             for Autonomous Diving and Steering of Unmanned Underwater Vehicles,
%             IEEE Journal of Ocean Engineering, JOE-18(3):327-339
%
% Author:    Trygve Lauvdal
% Date:      15-May-1994
% Revisions: 21-Jul-2001, Thor I. Fossen: minor change of notation, including speed U in call
%            24-Mar-2003, Gianluca Antonelli/Thor I. Fossen: corrected the integrals for 
%                         Cy, Cz, Cn and Cm by introducing a new variable temp.
%            24-Mar-2003, Gianluca Antonelli/Thor I. Fossen: corrected an error in the restoring 
%                         force Z, that is (W-B)*cos(theta)*cos(phi).
%            13-Nov-2014, Anand Sundaresan, Mwn*w*n was corrected to Mwn*w*u
%            01-Aug-2019, David Hansch, corrected bugs in the cross-flow drag formulae

% Check of input and state dimensions
if (length(x) ~= 12),error('x-vector must have dimension 12 !');end
if (length(ui) ~= 6),error('u-vector must have dimension 6 !');end

% Dimensional states
u   = x(1);  v     = x(2);  w   = x(3);
p   = x(4);  q     = x(5);  r   = x(6);
phi = x(10); theta = x(11); psi = x(12);

U = sqrt(u^2+v^2+w^2); % speed

% Rudder and propeller
max_ui(1)  = 20*pi/180;        % max value delta_r   (rad)
max_ui(2)  = 20*pi/180;        % max value delta_s   (rad)
max_ui(3)  = 20*pi/180;        % max value delta_b   (rad)
max_ui(4)  = 20*pi/180;        % max value delta_bp  (rad)
max_ui(5)  = 20*pi/180;        % max value delta_bs  (rad)
max_ui(6)  = 1500;             % max value n         (rpm)

% Parameters, hydrodynamic derivatives and main dimensions
c1 = cos(phi);
c2 = cos(theta);
c3 = cos(psi);
s1 = sin(phi);
s2 = sin(theta);
s3 = sin(psi);
t2 = tan(theta);

L   = 5.3;    g   = 9.8;
xG  = 0;      yG  = 0;      zG = 0.061;
xB  = 0;      yB  = 0;      zB = 0;
rho = 1000;   m   = 5454.54/(rho/2*L^3);
W   = 53400;  B   = 53400;
Ix  = 2038;   Iy  = 13587;  Iz  = 13587;
Ixy = -13.58; Iyz = -13.58; Ixz = -13.58;
Cdy = 0.5;    Cdz = 0.6;
Cy  = 0;      Cz  = 0;
Cm  = 0;      Cn  = 0;

r2 = rho*L^2/2;
r3 = rho*L^3/2;
r4 = rho*L^4/2;
r5 = rho*L^5/2;

Xpp   =  7.0e-3; Xqq    = -1.5e-2; Xrr   =  4.0e-3; Xpr   =  7.5e-4;
Xudot = -7.6e-3; Xwq    = -2.0e-1; Xvp   = -3.0e-3; Xvr   =  2.0e-2;
Xqds  =  2.5e-2; Xqdb2  = -1.3e-3; Xrdr  = -1.0e-3; Xvv   =  5.3e-2;
Xww   =  1.7e-1; Xvdr   =  1.7e-3; Xwds  =  4.6e-2; Xwdb2 =  0.5e-2;
Xdsds = -1.0e-2; Xdbdb2 = -4.0e-3; Xdrdr = -1.0e-2; Xqdsn =  2.0e-3;
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

% Rudder and shaft saturations
for i=1:1:6
   if abs(ui(i))>max_ui(i),ui(i)=sign(ui(i))*max_ui(i);end
end
 
% Control input (rudder and propeller)
delta_r  = ui(1);
delta_s  = ui(2);
delta_b  = ui(3);
delta_bp = ui(4);
delta_bs = ui(5);
n        = ui(6)/60*2*pi;

Cd0   = 0.00385;
prop  = 0.012*n/u; 
Xprop = Cd0*(abs(prop)*prop - 1);
Ct    = 0.008*L^2*abs(prop)*prop/2;
Ct1   = 0.008*L^2/2;
epsi  = -1 + sign(n)/sign(u)*(sqrt(Ct+1)-1)/(sqrt(Ct1+1)-1);

tau1 = r3*(Xrdr*u*r*delta_r + (Xqds*delta_s + Xqdb2*delta_bp +...
           Xqdb2*delta_bs)*u*q) +...
       r2*(Xvdr*u*v*delta_r + (Xwds*delta_s + Xwdb2*delta_bs + ...
           Xwdb2*delta_bp)*u*w + (Xdsds*delta_s^2 + ...
           Xdbdb2*delta_b^2 + Xdrdr*delta_r^2)*u^2) +...
       r3*Xqdsn*u*q*delta_s*epsi + r2*(Xwdsn*u*w*delta_s +...
       Xdsdsn*u^2*delta_s^2)*epsi + r2*u^2*Xprop; 
tau2 = r2*Ydr*u^2*delta_r;
tau3 = r2*u^2*(Zds*delta_s+Zdb2*delta_bs+Zdb2*delta_bp) + ...
       r3*Zqn*u*q*epsi+r2*(Zwn*u*w+Zdsn*u^2*delta_s)*epsi;
tau4 = r4*Kpn*u*p*epsi+r3*u^3*Kprop + ...
       r3*u^2*(Kdb2*delta_bp+Kdb2*delta_bs);
tau5 = r4*Mqn*u*q*epsi+r3*(Mwn*w*u+Mdsn*u^2*delta_s)*epsi+...
       r3*u^2*(Mds*delta_s+Mdb2*delta_bp+Mdb2*delta_bs);
tau6 = r3*u^2*Nprop + r3*u^2*Ndr*delta_r;

% Drag forces and moments assuming block shaped body
dxL = L/10;
xL  = 0; 
Ucf = sqrt((v+xL*r)^2+(w-xL*q)^2);

if ~(Ucf == 0)
    for xL = -L/2:dxL:L/2
        Ucf = sqrt((v+xL*r)^2+(w-xL*q)^2);
        temp = (0.5*0.6*(v+xL*r)^2+0.6*(w-xL*q)^2)*(v+xL*r)/Ucf; 
        Cy = Cy + dxL*temp;
    end
    
    for xL = -L/2:dxL:L/2
        Ucf = sqrt((v+xL*r)^2+(w-xL*q)^2);
        temp = (0.5*0.6*(v+xL*r)^2+0.6*(w-xL*q)^2)*(w-xL*q)/Ucf;
        Cz = Cz + dxL*temp;
    end
    
    for xL = -L/2:dxL:L/2
        Ucf = sqrt((v+xL*r)^2+(w-xL*q)^2);
        temp = (0.5*0.6*(v+xL*r)^2+0.6*(w-xL*q)^2)*(w+xL*q)/Ucf*xL;
        Cm = Cm + dxL*temp;
    end
    
    for xL = -L/2:dxL:L/2
        Ucf = sqrt((v+xL*r)^2+(w-xL*q)^2);
        temp = (0.5*0.6*(v+xL*r)^2+0.6*(w-xL*q)^2)*(v+xL*r)/Ucf*xL;
        Cn = Cn + dxL*temp;
    end
end

% Hydrodynamic forces and moments
X = r3*((m+Xvr)*v*r + (Xwq-m)*w*q + Xvp*v*p) + ... 
    r4*((m*xG/L+Xqq)*q^2 + (m*xG/L+Xrr)*r^2 - m*yG/L*p*q + ...
       (Xpr-m*zG/L)*p*r + Xpp*p^2 ) + ...
    r2*(Xvv*v^2 + Xww*w^2 )- (W - B)*sin(theta) + tau1;

Y = r2*(Yv*u*v + Yvw*v*w) + ...
    r3*(Yp*u*p + Yr*u*r +Yvq*v*q + Ywp*w*p + Ywr*w*r) +...
    r4*(Ypq*p*q + Yqr*q*r) + (W-B)*cos(theta)*sin(phi) - ...
    m*(rho/2*L^3)*(u*r -w*p + xG*p*q - yG*(p^2+r^2) + zG*q*r)+tau2-rho/2*Cy;

Z = r2*(Zw*w*u + Zvv*v^2) +...
    r3*(Zq*u*q + Zvp*v*p + Zvr*v*r) +...
    r4*(Zpp*p^2 + Zpr*p*r + Zrr*r^2) + ...
    (W-B)*cos(theta)*cos(phi) - ...
    m*(rho/2*L^3)*(v*p - u*q + xG*p*r + yG*q*r -zG*(p^2+q^2))+tau3+ rho/2*Cz;

K = r3*(Kv*u*v + Kvw*v*w) +...
    r4*(Kp*u*p + Kr*u*r + Kvq*v*q + Kwp*w*p + Kwr*w*r) +...
    r5*(Kpq*p*q + Kqr*q*r) +...
    (Iy-Iz)*q*r - Ixy*p*r - (r^2-q^2)*Iyz + Ixz*p*q - ...
    m*(rho/2*L^3)*(yG*(v*p-u*q) - zG*(u*r - w*p)) + ...
    (yG*W-yB*B)*c1*c2 - (zG*W-zB*B)*c2*s1 + tau4;

M = r3*(Muw*u*w + Mvv*v^2) +...
    r4*(Muq*u*q + Mvp*v*p + Mvr*v*r) +...
    r5*(Mpp*p^2 + Mpr*p*r + Mrr*r^2) - ...
    (xG*W-xB*B)*c1*c2 - (zG*W-zB*B)*s2 + ...
    (Iz-Ix)*p*r + Ixy*q*r - Iyz*p*q - (p^2-r^2)*Ixz + ...
    m*(rho/2*L^3)*(xG*(v*p-u*q) - zG*(w*q-v*r)) + tau5 - rho/2*Cm;

N = r3*(Nv*u*v + Nvw*v*w)+...
    r4*(Np*u*p + Nr*u*r + Nvq*v*q + Nwp*w*p + Nwr*w*r)+...
    r5*(Npq*p*q + Nqr*q*r) +...
    (xG*W-xB*B)*s1*c2 + (yG*W-yB*B)*s2 + ...
    (Ix-Iy)*p*q + (p^2-q^2)*Ixy + Iyz*p*r - Ixz*q*r -...
    m*(rho/2*L^3)*(xG*(u*r-w*p) - yG*(w*q-v*r)) + tau6 - rho/2*Cn ;

% Dimensional state derivatives ( xdot = in(M)*f(x) is expanded to avoid inv(M) on-line )
xdot = [1.662e-4*X+1.846e-10*Y+1.303e-7*Z+3.726e-9*K-1.132e-6*M+7.320e-10*N
        1.846e-10*X+1.052e-4*Y+3.843e-10*Z+9.638e-6*K-3.340e-9*M+2.368e-6*N
         1.303e-7*X+3.843e-10*Y+4.315e-5*Z+7.757e-9*K-2.357e-6*M+1.524e-9*N
          3.726e-9*X+9.638e-6*Y+7.757e-9*Z+2.431e-4*K-6.742e-8*M-7.740e-7*N
         -1.132e-6*X-3.340e-9*Y-2.357e-6*Z-6.742e-8*K+2.049e-5*M-1.324e-8*N
         7.320e-10*X+2.368e-6*Y+1.524e-9*Z-7.740e-7*K-1.324e-8*M+4.838e-5*N    
            c3*c2*u + (c3*s2*s1-s3*c1)*v + (s3*s1+c3*c1*s2)*w 
            s3*c2*u + (c1*c3+s1*s2*s3)*v + (c1*s2*s3-c3*s1)*w 
                     -s2*u + c2*s1*v + c1*c2*w 
                        p + s1*t2*q + c1*t2*r 
                            (c1*q - s1*r) 
                           s1/c2*q + c1/c2*r                               ];