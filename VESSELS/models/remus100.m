function [xdot,U,M] = remus100(x,ui,Vc,betaVc,w_c)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% The length of the Remus 100 AUV is 1.6 m, the cylinder diameter is 19 cm  
% and the mass of the vehicle is 31.9 kg. The maximum speed of 2.5 m/s is 
% obtained when the propeller runs at 1525 rpm in zero currents. The
% function returns the time derivative xdot of the state vector: 
%
%   x = [ u v w p q r x y z phi theta psi ]', alternatively 
%   x = [ u v w p q r x y z eta eps1 eps2 eps3 ]' 
%
% in addition to the speed U in m/s (optionally). The state vector can be 
% of dimension 12 (Euler angles) or 13 (unit quaternions):
%
%   u:       Surge velocity          (m/s)
%   v:       Sway velocity           (m/s)
%   w:       Heave velocity          (m/s)
%   p:       Roll rate               (rad/s)
%   q:       Pitch rate              (rad/s)
%   r:       Yaw rate                (rad/s)
%   x:       North position          (m)
%   y:       East position           (m)
%   z:       Downwards position      (m)
%   phi:     Roll angle              (rad)       
%   theta:   Pitch angle             (rad)
%   psi:     Yaw angle               (rad)
% 
% For the unit quaternion representation, the last three arguments of the 
% x-vector, the Euler angles (phi, theta, psi), are replaced by the unit 
% quaternion quat = [eta, eps1, eps2, eps3]'. This increases the dimension of
%  the state vector from 12 to 13.
%
% The control inputs are one tail rudder, two stern planes and a single-screw 
% propeller:
%
%   ui = [ delta_r delta_s n ]'  where
%
%    delta_r:   Rudder angle (rad)
%    delta_s:   Stern plane angle (rad) 
%    n:         Propeller revolution (RPM)
%
% The arguments Vc (m/s), betaVc (rad), w_c (m/s) are optional arguments for 
% ocean currents
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), w_c ]  
% 
% Example usage: 
%   [~,~,M] = remus100()                            : Return the 6x6 mass matrix M
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc,w_c) : 3-D ocean currents
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc)     : 2-D ocean currents
%   [xdot,U] = remus100(x,ui)                       : No ocean currents
%   xdot = remus100(x,ui)                           : No ocean currents
%
% Author:    Thor I. Fossen
% Date:      2021-05-27
% Revisions:
%   2021-08-24  Ocean currents are now expressed in NED
%   2021-10-21  imlay61.m is called using the relative velocity
%   2021-12-30  Added the time derivative of the current velocity
%   2022-02-01  Updated lift and drag forces
%   2022-05-06  Calibration of drag and propulsion forces using data from 
%               Allen et al. (2000)
%   2022-06-08  Added compatibility for unit quaternions in addition to the 
%               Euler angle representation
%   2022-10-16  Added vertical currents
%   2023-05-02  Corrected the rudder area A_r
%   2023-10-07  Scaled down the propeller roll-induced moment
%   2024-02-09  Updated rudder and stern-plane areas
%   2024-02-13  Calibration of the model parameters
%   2024-06-23  Corrected sign of tau(5). A positive delta_s will result in a
%               negative pitch and negative Z_s. Hence, tau(5) = -x_s * Z_s 
%               when x_s < 0.
%   2025-04-25 Added empty call: [~,~,M] = remus100(), and minor bug fixes.
%
% References: 
%   B. Allen, W. S. Vorus and T. Prestero, "Propulsion system 
%       performance enhancements on REMUS AUVs," OCEANS 2000 MTS/IEEE 
%       Conference and Exhibition. Conference Proceedings, 2000, 
%       pp. 1869-1873 vol.3, doi: 10.1109/OCEANS.2000.882209.
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 2nd. Edition, Wiley. URL: www.fossen.biz/wiley   

if nargin == 0  % [~,~,M] = remus100()
    x = zeros(12,1); ui = zeros(3,1); Vc = 0; betaVc = 0; w_c = 0;
end

if (nargin == 2), Vc = 0; betaVc = 0; w_c = 0; end % No ocean currents
if (nargin == 4), w_c = 0; end % No vertical ocean currents

if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end
if (length(x) ~= 12 && length(x) ~= 13)
    error('x-vector must have dimension 12 or 13'); 
end

% Constants
mu = 63.446827;         % Lattitude for Trondheim, Norway (deg)
g_mu = gravity(mu);     % Gravity vector (m/s2)
rho = 1026;             % Density of water (m/s2)

% 6x1 velocity vector, yaw angle and control inputs
nu = x(1:6); 
if length(x) == 12 % Euler angles
    psi = x(12);
else % Convert unit quaternion to yaw angle
    psi = atan2(2*(x(10)*x(13) + x(11)*x(12)), 1 - 2*(x(12)^2 + x(13)^2));
end

% Amplitude saturation of the control signals
delta_max = deg2rad(20); % Maximum rudder and stern angles (rad)
n_max = 1525;            % Maximum propeller speed (RPM)

% Amplitude saturation of the control signals
delta_r = sat(ui(1), delta_max);    % Saturated tail rudder (rad)
delta_s = sat(ui(2), delta_max);    % Saturated Stern plane (rad)
n = sat(ui(3),n_max) / 60;          % Saturated propeller speed (rps)

% Ocean currents expressed in BODY
u_c = Vc * cos( betaVc - psi );                               
v_c = Vc * sin( betaVc - psi );   

nu_c = [u_c v_c w_c 0 0 0]'; % Ocean current velocities
Dnu_c = [nu(6)*v_c -nu(6)*u_c 0 0 0 0]'; % Time derivative of nu_c

% Relative velocities/speed, angle of attack and vehicle speed
nu_r = nu - nu_c;                                 % Relative velocity
alpha = atan2( nu_r(3), nu_r(1) );                % Angle of attack (rad)
U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );  % Relative speed (m/s)
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );         % Apeed (m/s)

% AUV model parameters; Fossen (2021, Section 8.4.2) and Allen et al. (2000)
L_auv = 1.6;             % AUV length (m)
D_auv = 0.19;            % AUV diamater (m)
S = 0.7 * L_auv * D_auv; % Planform area S = 70% of rectangle L_auv * D_auv
a = 1.0096 * L_auv/2;    % Scaled spheroid semi-axes a and b to obtain m = 31.9 kg
b = 1.0096 * D_auv/2;                   
r44 = 0.3;               % Added moment of inertia in roll: A44 = r44 * Ix
r_bG = [ 0 0 0.02 ]';    % CG w.r.t. to the CO
r_bB = [ 0 0 0 ]';       % CB w.r.t. to the CO

% Parasitic drag coefficient CD_0, i.e. zero lift and alpha = 0
% F_drag = 0.5 * rho * Cd * (pi * b^2)   
% F_drag = 0.5 * rho * CD_0 * S
Cd = 0.42;                              % From Allen et al. (2000)
CD_0 = Cd * pi * b^2 / S;

% Propeller coeffs. KT and KQ are computed as a function of advance no.
% Ja = Va/(n*D_prop) where Va = (1-w)*U = 0.944 * U; Allen et al. (2000)
D_prop = 0.14;   % Propeller diameter corresponding to 5.5 inches
t_prop = 0.1;    % Thrust deduction number
Va = 0.944 * U;  % Advance speed (m/s)

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
if n > 0   
    % Forward thrust
    X_prop = rho * D_prop^4 * (... 
        KT_0 * abs(n) * n + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n) );        
    K_prop = rho * D_prop^5 * (...
        KQ_0 * abs(n) * n + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n) );                 
else      
    % Reverse thrust (braking)   
    X_prop = rho * D_prop^4 * KT_0 * abs(n) * n; 
    K_prop = rho * D_prop^5 * KQ_0 * abs(n) * n;            
end            

S_fin = 0.00665;         % Fin area

% Tail rudder
CL_delta_r = 0.5;        % Rudder lift coefficient (-)
A_r = 2 * S_fin;         % Rudder area (m2)
x_r = -a;                % Rudder x-position (m)

% Stern plane (double)
CL_delta_s = 0.7;        % Stern-plane lift coefficient (-)
A_s = 2 * S_fin;         % Stern-plane area (m2)
x_s = -a;                % Stern-plane z-position (m)

% Low-speed linear damping matrix parameters
T1 = 20;                 % Time constant in surge (s)
T2 = 20;                 % Time constant in sway (s)
zeta4 = 0.3;             % Relative damping ratio in roll
zeta5 = 0.8;             % Relative damping ratio in pitch
T6 = 1;                  % Time constant in yaw (s)

% Rigid-body mass and hydrodynamic added mass
[MRB,CRB] = spheroid(a,b,nu(4:6),r_bG);
[MA,CA] = imlay61(a, b, nu_r, r44);

% CA-terms in roll, pitch and yaw can destabilize the model if quadratic
% rotational damping is missing. These terms are assumed to be zero
CA(5,3) = 0; CA(3,5) = 0;  % Quadratic velocity terms due to pitching
CA(5,1) = 0; CA(1,5) = 0;  
CA(6,1) = 0; CA(1,6) = 0;  % Munk moment in yaw 
CA(6,2) = 0; CA(2,6) = 0;

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1); W = m * g_mu; B = W;

% Dissipative forces and moments
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);
D(1,1) = D(1,1) * exp(-3*U_r);   % Vanish at high speed where quadratic
D(2,2) = D(2,2) * exp(-3*U_r);   % Drag and lift forces dominates

tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r);

% Kinematics
if (length(x) == 13)
    [J,R] = quatern(x(10:13));
else
    [J,R] = eulerang(x(10),x(11),x(12));
end

% Restoring forces and moments
g = gRvect(W,B,R,r_bG,r_bB);

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
tau(4) = K_prop / 10; % Scaled down by a factor of 10 to match exp. results
tau(5) = -x_s * Z_s;  
tau(6) = x_r * Y_r;

% State-space model
xdot = [ Dnu_c + M \ ...
            (tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)
         J * nu ]; 