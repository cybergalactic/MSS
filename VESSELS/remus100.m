function [xdot,U] = remus100(x,ui,v_current)
% [xdot,U] = remus100(x,ui,v_current) returns the time derivative of the 
% state vector: x = [ u v w p q r x y z phi theta psi ]' and speed U in m/s  
% (optionally) for the Remus 100 autonomous underwater vehicle (AUV). The 
% length of the AUV is L = 1.7 m, while the state vector is defined as:
%
%  u:       surge velocity          (m/s)
%  v:       sway velocity           (m/s)
%  w:       heave velocity          (m/s)
%  p:       roll rate               (rad/s)
%  q:       pitch rate              (rad/s)
%  r:       yaw rate                (rad/s)
%  x:       North position          (m)
%  y:       East position           (m)
%  z:       downwards position      (m)
%  phi:     roll angle              (rad)
%  theta:   pitch angle             (rad)
%  psi:     yaw angle               (rad)
%
% The control inputs are:
%
%  ui = [ delta_r delta_s n ]'  where
%
%    delta_r:   rudder angle (rad)
%    delta_s:   aft stern plane (rad) 
%    n:         propeller revolution (rpm)
%
% The last argument v_current is an optional argument for ocean current 
% velocities v_current = [u_c v_c w_c]' expressed in the BODY frame.
%
% Author:    Thor I. Fossen
% Date:      6 May 2021
% Revisions: THIS IS A TEMPORARY MODEL - USE WITH CARE

% Check of input and state dimensions
if (length(x) ~= 12),error('x-vector must have dimension 12!'); end
if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end

% Constants
mu = 63.446827;         % Lattitude for Trondheim (deg)
g_mu = gravity(mu);     % gravity vector (m/s2)
rho = 1026;             % density of water (m/s2)

% State vectors and control inputs
nu = x(1:6); 
eta = x(7:12);
delta_r = ui(1);        % tail rudder (rad)
delta_s = ui(2);        % stern plane (rad)
n = ui(3)/60;           % propeller revolution (rps)

% Amplitude saturation of control signals
max_ui = [20*pi/180 20*pi/180  1500/60]';  % deg, deg, rps

% Relative velocities, speed and angle of attack
if (nargin == 2), v_current = [0 0 0]'; end
nu_r = nu - [v_current' 0 0 0]';                  % relative velocity vevctor
alpha = atan2( nu_r(3), nu_r(1) );                % angle of attack (rad)
U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );  % relative speed (m/s)
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );         % speed (m/s)

% AUV model parameters (Section 8.4.2)
L_auv = 1.7;                        % AUV length (m)
D_auv = 0.19;                       % AUV diamater (m)
CD_0 = 0.01;                        % parasitic drag
S = 0.7 * L_auv * D_auv;            % S = 70% of rectangle L_auv * D_auv
a = L_auv/2;                        % semi-axes
b = D_auv/2;                  
r_bg = [ 0 0 0.025 ]';               % CG w.r.t. to the CO
r_bb = [ 0 0 0 ]';                  % CB w.r.t. to the CO

% Added moment of inertia in roll
r44 = 0.3;               % A44 = r44 * Ix

% Propeller data
D_prop = 0.055;          % propeller diameter
KT = 0.4739;             % [KT, KQ] = wageningen(0,1,0.8,3)
KQ = 0.0730;

% Tail rudder
CL_delta_r = 0.5;        % rudder lift coefficient
A_r = 0.10 * 0.05;       % rudder area (m2)
x_r = -a;                % rudder x-position (m)

% Stern plane
CL_delta_s = 0.7;        % stern-plane lift coefficient
A_s = 0.10 * 0.05;       % stern-plane area (m2)
x_s = -a;                % stern-plane z-position (m)

% Low-speed linear damping matrix parameters: D * exp(-3 * U_r)
T1 = 20;                 % time constant in surge (s)
T2 = 20;                 % time constant in sway (s)
zeta4 = 0.3;             % relative damping ratio in roll
zeta5 = 0.8;             % relative damping ratio in pitch
T6 = 5;                  % time constant in yaw (s)

% mass and added mass
[MRB,CRB] = spheroid(a,b,nu(4:6),r_bg);
[MA,CA] = imlay61(a, b, nu, r44);

% nonlinear quadratic velocity terms in pitch and yaw (Munk moments) 
% are cancelled since only linear damping is used
CA(5,1) = 0;   
CA(5,4) = 0;
CA(6,1) = 0;
CA(6,2) = 0;

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1); W = m*g_mu; B = W;

% dissipative forces and moments
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bg' r_bb']);
D(1,1) = D(1,1) * exp(-3*U_r);   % vanish at high speed where crossflow
D(2,2) = D(2,2) * exp(-3*U_r);   % drag and lift/drag dominates
D(6,6) = D(6,6) * exp(-3*U_r);

tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r);

% restoring forces and moments
g = gvect(W,B,eta(5),eta(4),r_bg,r_bb);

% kinematics
J = eulerang(eta(4),eta(5),eta(6));

% amplitude saturation of the control signals
if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end

% control forces and moments
X_prop = rho * D_prop^4 * KT * abs(n) * n;  % propeller thrust 
K_prop = rho * D_prop^5 * KQ * abs(n) * n;  % propeller-induced roll moment

X_r = -0.5 * rho * U_r^2 * A_r * CL_delta_r * delta_r^2;  % rudder drag
Y_r = -0.5 * rho * U_r^2 * A_r * CL_delta_r * delta_r;    % rudder sway force
N_r = x_r * Y_r;                                          % rudder yaw moment

X_s = -0.5 * rho * U_r^2 * A_s * CL_delta_s * delta_s^2;  % stern-plane drag
Z_s = -0.5 * rho * U_r^2 * A_s * CL_delta_s * delta_s;    % stern-plane heave force
M_s = -x_s * Z_s;                                         % stern-plane pitch moment

tau = zeros(6,1);                                % generalized force vector
tau(1) = X_prop + X_r + X_s;
tau(2) = Y_r;
tau(3) = Z_s;
tau(4) = K_prop;
tau(5) = M_s;
tau(6) = N_r;

% state-space model
xdot = [...
    M  \ (tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)
    J * nu ];

