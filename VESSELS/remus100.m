function [xdot,U] = remus100(x,ui,Vc,betaVc,w_c)
% The length of the Remus 100 AUV is 1.6 m, the cylinder diameter is 19 cm  
% and the mass of the vehicle is 31.9 kg. The maximum speed of 2.5 m/s is 
% obtained when the propeller runs at 1525 rpm in zero currents. The
% function calls are:
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc,w_c)  3-D ocean currents
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc)      horizontal ocean currents
%   [xdot,U] = remus100(x,ui)                        no ocean currents
% The function returns the time derivative xdot of the state vector: 
%   x = [ u v w p q r x y z phi theta psi ]',     alternatively 
%   x = [ u v w p q r x y z eta eps1 eps2 eps3 ]' 
% in addition to the speed U in m/s (optionally). The state vector can be 
% of dimension 12 (Euler angles) or 13 (unit quaternions):
%
%   u:       surge velocity          (m/s)
%   v:       sway velocity           (m/s)
%   w:       heave velocity          (m/s)
%   p:       roll rate               (rad/s)
%   q:       pitch rate              (rad/s)
%   r:       yaw rate                (rad/s)
%   x:       North position          (m)
%   y:       East position           (m)
%   z:       downwards position      (m)
%   phi:     roll angle              (rad)       
%   theta:   pitch angle             (rad)
%   psi:     yaw angle               (rad)
% 
% For the unit quaternion representation, the last three arguments of the 
% x-vector, the Euler angles (phi, theta, psi), are replaced by the unit 
% quaternion q = [eta, eps1, eps2, eps3]'. This increases the dimension of 
% the state vector from 12 to 13.
%
% The control inputs are one tail rudder, two stern planes and a single-screw 
% propeller:
%
%   ui = [ delta_r delta_s n ]'  where
%
%    delta_r:   rudder angle (rad)
%    delta_s:   stern plane angle (rad) 
%    n:         propeller revolution (rpm)
%
% The arguments Vc (m/s), betaVc (rad), w_c (m/s) are optional arguments for 
% ocean currents
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), w_c ]  
% 
% Author:    Thor I. Fossen
% Date:      27 May 2021
% Revisions: 24 Aug 2021  Ocean currents are now expressed in NED 
%            21 Oct 2021  imlay61.m is called using the relative velocity
%            30 Dec 2021  Added the time derivative of the current velocity
%            01 Feb 2022  Updated lift and drag forces
%            06 May 2022  Calibration of drag and propulsion forces using
%                         data from Allen et al. (2000)
%            08 May 2022  Added compability for unit quaternions in 
%                         addition to the Euler angle representation
%            16 Oct 2022  Added vertical currents
%
% Refs: 
%      B. Allen, W. S. Vorus and T. Prestero, "Propulsion system 
%           performance enhancements on REMUS AUVs," OCEANS 2000 MTS/IEEE 
%           Conference and Exhibition. Conference Proceedings, 2000, 
%           pp. 1869-1873 vol.3, doi: 10.1109/OCEANS.2000.882209.
%      T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%           Motion Control. 2nd. Edition, Wiley. URL: www.fossen.biz/wiley   

if (nargin == 2), Vc = 0; betaVc = 0; w_c = 0; end  % no ocean currents
if (nargin == 4), w_c = 0; end             % no vertical ocean currents

if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end
if (length(x) ~= 12 && length(x) ~= 13)
    error('x-vector must have dimension 12 or 13'); 
end

% Constants
mu = 63.446827;         % Lattitude for Trondheim, Norway (deg)
g_mu = gravity(mu);     % gravity vector (m/s2)
rho = 1026;             % density of water (m/s2)

% State vectors and control inputs
nu = x(1:6); 
eta = x(7:12);
delta_r = ui(1);        % tail rudder (rad)
delta_s = ui(2);        % stern plane (rad)
n = ui(3)/60;           % propeller revolution (rps)

% Ocean currents expressed in BODY
u_c = Vc * cos( betaVc - eta(6) );                               
v_c = Vc * sin( betaVc - eta(6) );   

nu_c = [u_c v_c w_c 0 0 0]';                  % ocean current velocities
Dnu_c = [nu(6)*v_c -nu(6)*u_c 0 0 0 0]';    % time derivative of nu_c

% Amplitude saturation of rudder angle, stern plane and propeller revolution
n_max = 1525;                                % maximum propeller rpm
max_ui = [30*pi/180 30*pi/180  n_max/60]';   % deg, deg, rps

% Relative velocities/speed, angle of attack and vehicle speed
nu_r = nu - nu_c;                                 % relative velocity
alpha = atan2( nu_r(3), nu_r(1) );                % angle of attack (rad)
U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );  % relative speed (m/s)
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );         % speed (m/s)

% AUV model parameters; Fossen (2021, Section 8.4.2) and Allen et al. (2000)
L_auv = 1.6;             % AUV length (m)
D_auv = 0.19;            % AUV diamater (m)
S = 0.7 * L_auv * D_auv; % planform area S = 70% of rectangle L_auv * D_auv
a = L_auv/2;             % spheroid semi-axes a and b
b = D_auv/2;                  
r44 = 0.3;               % added moment of inertia in roll: A44 = r44 * Ix
r_bg = [ 0 0 0.02 ]';    % CG w.r.t. to the CO
r_bb = [ 0 0 0 ]';       % CB w.r.t. to the CO

% Parasitic drag coefficient CD_0, i.e. zero lift and alpha = 0
% F_drag = 0.5 * rho * Cd * (pi * b^2)   
% F_drag = 0.5 * rho * CD_0 * S
Cd = 0.42;                              % from Allen et al. (2000)
CD_0 = Cd * pi * b^2 / S;

% Propeller coeffs. KT and KQ are computed as a function of advance no.
% Ja = Va/(n*D_prop) where Va = (1-w)*U = 0.944 * U; Allen et al. (2000)
D_prop = 0.14;   % propeller diameter corresponding to 5.5 inches
t_prop = 0.1;    % thrust deduction number
Va = 0.944 * U;  % advance speed (m/s)

% Ja_max = 0.944 * 2.5 / (0.14 * 1525/60) = 0.6632
Ja_max = 0.6632;
        
% Single-screw propeller with 3 blades and blade-area ratio = 0.718.    
% >> [KT_0, KQ_0] = wageningen(0,1,0.718,3)
KT_0 = 0.4566;
KQ_0 = 0.0700;
% >> [KT_max, KQ_max] = wageningen(0.6632,1,0.718,3) 
KT_max = 0.1798;
KQ_max = 0.0312;
        
% Propeller thrust and propeller-induced roll moment
% Linear approximations for positive Ja values
% KT ~= KT_0 + (KT_max-KT_0)/Ja_max * Ja   
% KQ ~= KQ_0 + (KQ_max-KQ_0)/Ja_max * Ja  
      
if n > 0   % forward thrust
    X_prop = rho * D_prop^4 * (... 
        KT_0 * abs(n) * n + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n) );        
    K_prop = rho * D_prop^5 * (...
        KQ_0 * abs(n) * n + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n) );           
            
else    % reverse thrust (braking)
        
    X_prop = rho * D_prop^4 * KT_0 * abs(n) * n; 
    K_prop = rho * D_prop^5 * KQ_0 * abs(n) * n;
            
end            

% Tail rudder (single)
CL_delta_r = 0.5;        % rudder lift coefficient (-)
A_r = 0.10 * 0.05;       % rudder area (m2)
x_r = -a;                % rudder x-position (m)

% Stern plane (double)
CL_delta_s = 0.7;        % stern-plane lift coefficient (-)
A_s = 2 * 0.10 * 0.05;   % stern-plane area (m2)
x_s = -a;                % stern-plane z-position (m)

% Low-speed linear damping matrix parameters
T1 = 20;                 % time constant in surge (s)
T2 = 20;                 % time constant in sway (s)
zeta4 = 0.3;             % relative damping ratio in roll
zeta5 = 0.8;             % relative damping ratio in pitch
T6 = 5;                  % time constant in yaw (s)

% Rigid-body mass and hydrodynamic added mass
[MRB,CRB] = spheroid(a,b,nu(4:6),r_bg);
[MA,CA] = imlay61(a, b, nu_r, r44);

% Nonlinear quadratic velocity terms in pitch and yaw. Munk moments 
% are set to zero since only linear rotational damping is used in the model
CA(5,1) = 0;   
CA(5,4) = 0;
CA(6,1) = 0;
CA(6,2) = 0;

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1); W = m*g_mu; B = W;

% Dissipative forces and moments
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bg' r_bb']);
D(1,1) = D(1,1) * exp(-3*U_r);   % vanish at high speed where quadratic
D(2,2) = D(2,2) * exp(-3*U_r);   % drag and lift forces dominates
D(6,6) = D(6,6) * exp(-3*U_r);

tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r);

% Kinematics
if (length(x) == 13)
    [J,R] = quatern(x(10:13));
else
    [J,R] = eulerang(x(10),x(11),x(12));
end

% Restoring forces and moments
g = gRvect(W,B,R,r_bg,r_bb);

% Amplitude saturation of the control signals
if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end

% Horizontal- and vertical-plane relative speed
U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );  
U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );  

% Rudder and stern-plane drag
X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2; 
X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;

% Rudder sway force 
Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;

% Stern-plane heave force
Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;

% Generalized propulsion force vector
tau = zeros(6,1);                                
tau(1) = (1-t_prop) * X_prop + X_r + X_s;
tau(2) = Y_r;
tau(3) = Z_s;
tau(4) = K_prop;
tau(5) = x_s * Z_s;
tau(6) = x_r * Y_r;

% State-space model
xdot = [ Dnu_c + M \ ...
            (tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)
         J * nu ]; 