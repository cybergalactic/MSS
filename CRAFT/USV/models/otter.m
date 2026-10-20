function [xdot,U,M,B_prop,n_min,n_max] = otter(x,n,mP,r_bP,V_c,beta_c)
% [xdot,U,M,B_prop,n_min,n_max] = otter(x,n,mP,r_bP,V_c,beta_c) returns the
% time derivative xdot of the state vector, the speed U in m/s, the 6x6
% mass matrix M, the 2x2 propeller input matrix B_prop, and the physical
% propeller shaft-speed limits n_min and n_max.
%
%  x = [ u v w p q r x y z phi theta psi ]'
% 
% for the Maritime Robotics Otter USV, see www.maritimerobotics.com. 
% The length of the USV is L = 2.0 m, while the state vector is defined as:
%
%  u:     Surge velocity          (m/s)
%  v:     Sway velocity           (m/s)
%  w:     Heave velocity          (m/s)
%  p:     Roll velocity           (rad/s)
%  q:     Pitch velocity          (rad/s)
%  r:     Yaw velocity            (rad/s)
%  x:     Position in x direction (m)
%  y:     Position in y direction (m)
%  z:     Position in z direction (m)
%  phi:   Roll angle              (rad)
%  theta: Pitch angle             (rad)
%  psi:   Yaw angle               (rad)
%
% The other inputs are:
%
%  n = [ n(1) n(2) ]' where
%    n(1): Propeller shaft speed, left (rad/s)
%    n(2): Propeller shaft speed, right (rad/s)
%
%  mP: Payload mass (kg), maximum 45 kg and minimum 0 kg
%  r_bP = [xP, yP, zP]' (m) is the payload position from CO in body axes
%  V_c:     Ocean current speed (m/s)
%  beta_c:  Ocean current direction (rad)
%
% Position vectors r_bG = [xG, yG, zG]' and r_bP = [xP, yP, zP]' are
% measured from the body-frame origin CO to CG and P, respectively.
% The coordinate origin CO is chosen midships on the centerline. The 
% draft is T = 13.4 cm when the payload mP = 0. This corresponds to z = 0. 
% Choosing mP = 25 kg gives T = 19.5 cm and z = 0.04 m in steady state.
% Adding a non-zero payload mP will change the steady-state z value since 
% the payload is added as an external force g_0 in the code.
%
% Example usage:
%
%   [~,~,M,B_prop,n_min,n_max] = otter()
%       Return the 6x6 mass matrix M, the 2x2 input matrix B_prop, and the
%       propeller speed limits n_min and n_max.
%
%   [xdot,U] = otter(x,n,mP,r_bP)
%       Return xdot and U, no ocean currents.
%
%   xdot = otter(x,n,mP,r_bP)
%       Return xdot, no ocean currents.
%
%   xdot = otter(x,n,mP,r_bP,V_c,beta_c)
%       Return xdot, 2-D currents.
%
% M-file Simulators:
%   SIMotter.m : Script demonstrating 3-D ALOS path-following control.
%   exOtter.m  : Example script demonstrating the 5-state EKF for estimation 
%                of COG, COG, and course rate when using a course autopilot.
% 
% Simulink Simulators:
%   demoOtterUSVHeadingControl.slx : 3-D ALOS path-following control.
%
% Author:    Thor I. Fossen
% Date:      2019-07-17
% Revisions:
%   2021-04-25 : Added call to new function crossFlowDrag.m..
%   2021-07-22 : Added a new state for the trim moment.
%   2021-12-17 : New method Xudot = -addedMassSurge(m,L,rho).
%   2023-03-28 : Trim state is replaced by payload mP and r_bP.
%   2023-10-14 : Added ocean current acceleration terms.
%   2024-02-03 : Recalibration of damping terms.
%   2024-02-28 : Corrected the trim condition for eta.
%   2024-06-05 : Added two new output arguments for M and B_prop.
%   2026-04-20 : Correct payload lever arm moment of inertia.
%   2026-06-27 : Added output arguments for propeller speed limits n_min and n_max.
%   2026-10-07 : Corrected KB to use the pontoon waterplane area (E. Krizman)
%   2026-10-09 : Use loaded displacement mass in the surge added-mass estimate.
%   2026-10-14 : Calibrated damping using experimental time-series by blending 
%                linear and quadratic surge damping. For sway, T_sway is changed 
%                from 1.0 to 0.5 seconds to match experimental time responses. 

if nargin == 0
    x = zeros(12,1); n = zeros(2,1); mP=25; r_bP = zeros(3,1); V_c=0; beta_c=0;
end 

% Check of input and state dimensions
if (length(x) ~= 12),error('x vector must have dimension 12!'); end
if (length(n) ~= 2),error('n vector must have dimension 2!'); end

% Main data
constant = mssConstants();
g = constant.g;
rho = constant.rho_water;
L = 2.0;              % Length (m)
B = 1.08;             % Beam (m)
m = 55.0;             % Mass (kg)
r_bG = [0.2 0 -0.2]'; % Hull CG position from CO in body axes (m)
R44 = 0.4 * B;        % Radii of gyrations (m)
R55 = 0.25 * L;
R66 = 0.25 * L;
T_surge = 5;          % Time constant in surge (s)
T_sway = 0.5;         % Time constant in sway (s)
T_yaw = 1;            % Time constant in yaw (s)
Umax = 6 * 0.5144;    % 6 knots maximum forward speed (m/s)

% Data for one pontoon
B_pont  = 0.25;       % Beam of one pontoon (m)
y_pont  = 0.395;      % Distance from centerline to waterline area center (m)
Cw_pont = 0.75;       % Waterline area coefficient (-)
Cb_pont = 0.4;       % Block coefficient, computed from m = 55 kg

% State and current variables
nu = x(1:6);  nu1 = x(1:3); nu2 = x(4:6);   % Velocity vectors
eta = x(7:12);                              % Position vector
U = sqrt(nu(1)^2 + nu(2)^2 + nu(3)^2);      % Speed
u_c = V_c * cos(beta_c - eta(6));           % Current surge velocity
v_c = V_c * sin(beta_c - eta(6));           % Current sway velocity
nu_c = [u_c v_c 0 0 0 0 ]';                 % Current velocity vector
nu_r = nu - nu_c;                           % Relative velocity vector
nu_c_dot = [-Smtrx(nu2) * nu_c(1:3)         % Current acceleration vector
             zeros(3,1)              ];

% Inertia dyadic, volume displacement, and draft
nabla = (m+mP)/rho;                         % Volume
T = nabla / (2 * Cb_pont * B_pont * L);     % Draft
Ig_CG = m * diag([R44^2, R55^2, R66^2]);    % Only hull in the CG
r_bG_hull = r_bG;   
r_bG_total = (m*r_bG_hull + mP*r_bP)/(m+mP);
r_hull = r_bG_hull - r_bG_total;
r_payload = r_bP - r_bG_total;
Ig = Ig_CG - m * Smtrx(r_hull)^2 - mP * Smtrx(r_payload)^2;
r_bG = r_bG_total;

% Experimental propeller data including lever arms
l1 = -y_pont;                           % Lever arm, left propeller (m)
l2 = y_pont;                            % Lever arm, right propeller (m)
% The calibrated coefficients describe the combined twin-propeller thrust;
% division by two gives the per-propeller coefficients used in the thrust law.
k_pos = 0.02216/2;                      % Positive thrust, one propeller
k_neg = 0.01289/2;                      % Negative thrust, one propeller
n_max =  sqrt((0.5*24.4 * g)/k_pos);    % Maximum propeller rev. (rad/s)
n_min = -sqrt((0.5*13.6 * g)/k_neg);    % Minimum propeller rev. (rad/s)

% MRB and CRB (Fossen 2027)
I3 = eye(3);
O3 = zeros(3,3);

MRB_CG = [ (m+mP) * I3  O3
           O3           Ig ];
CRB_CG = [ (m+mP) * Smtrx(nu2)         O3
           O3               -Smtrx(Ig*nu2)  ];

H = Hmtrx(r_bG);              % Transform MRB and CRB from the CG to the CO 
MRB = H' * MRB_CG * H;
CRB = H' * CRB_CG * H;

% Hydrodynamic added mass (best practice)
Xudot = -addedMassSurge(m+mP,L,rho);
Yvdot = -1.5 * m;
Zwdot = -1.0 * m;
Kpdot = -0.2 * Ig(1,1);
Mqdot = -0.8 * Ig(2,2);
Nrdot = -1.7 * Ig(3,3);

MA = -diag([Xudot, Yvdot, Zwdot, Kpdot, Mqdot, Nrdot]);   
CA  = m2c(MA, nu_r);

% Uncomment to cancel the Munk moment in yaw, if stability problems
% CA(6,1) = 0; 
% CA(6,2) = 0; 
% CA(1,6) = 0;
% CA(2,6) = 0;

% System mass and Coriolis-centripetal matrices
M = MRB + MA;
C = CRB + CA;

% Hydrostatic quantities (Fossen 2027)
Aw_pont = Cw_pont * L * B_pont;    % Waterline area, one pontoon 
I_T = 2 * (1/12)*L*B_pont^3 * (6*Cw_pont^3/((1+Cw_pont)*(1+2*Cw_pont)))...
    + 2 * Aw_pont * y_pont^2;
I_L = 0.8 * 2 * (1/12) * B_pont * L^3;
KB = (1/3) * (5*T/2 - 0.5*nabla/Aw_pont);
BM_T = I_T/nabla;       % BM values
BM_L = I_L/nabla;
KM_T = KB + BM_T;       % KM values
KM_L = KB + BM_L;
KG = T - r_bG(3);
GM_T = KM_T - KG;       % GM values
GM_L = KM_L - KG;

G33 = rho * g * (2 * Aw_pont);      % Spring stiffness
G44 = rho * g * nabla * GM_T;
G55 = rho * g * nabla * GM_L;

G_CF = diag([0 0 G33 G44 G55 0]);   % Spring stiffness matrix in the CF
LCF = -0.2;
H = Hmtrx([LCF 0 0]);               % Transform G_CF from the CF to the CO 
G = H' * G_CF * H;

% Natural frequencies
w3 = sqrt( G33/M(3,3) );         
w4 = sqrt( G44/M(4,4) );
w5 = sqrt( G55/M(5,5) );

% Damping terms (hydrodynamic derivatives)
Xu = -M(1,1) / T_surge;       % Specified using the time constant in surge
Xuu = -24.4 * g / Umax^2;     % Specified using the maximum speed  
Yv = -M(2,2) / T_sway;        % Specified using the time constant in sway
Zw = -2 * 0.3 * w3 * M(3,3);  % Specified using relative damping factors
Kp = -2 * 0.2 * w4 * M(4,4);
Mq = -2 * 0.4 * w5 * M(5,5);
Nr = -M(6,6) / T_yaw;         % Specified using the time constant in T_yaw

% 2-DOF constant input matrix B_prop for the propellers in surge and yaw,
% accessible by: [~,~,M,B_prop,n_min,n_max] = otter()
if nargin == 0
    B_prop = k_pos * [...
        1 1                 
        y_pont -y_pont ];
else
    B_prop = [];
end

% Control forces and moments, with saturated propeller speed
n = satlim(n, n_min, n_max);

Thrust = zeros(2,1);
for i = 1:2
    if n(i) > 0
        Thrust(i) = k_pos * n(i) * abs(n(i)); % Positive thrust (N)
    else
        Thrust(i) = k_neg * n(i) * abs(n(i)); % Negative thrust (N)
    end
end

% Control forces and moments
tau = [Thrust(1) + Thrust(2) 0 0 0 0 -l1 * Thrust(1) - l2 * Thrust(2) ]';

% Damping using relative velocities + nonlinear surge and yaw damping
Xh = Xu  * exp(-3 * abs(nu_r(1))) * nu_r(1) + ...
    Xuu * abs(nu_r(1)) * nu_r(1); % Blending linear and quadratic Xu and Xuu
Yh = Yv * nu_r(2); 
Zh = Zw * nu_r(3);
Kh = Kp * nu_r(4);
Mh = Mq * nu_r(5);
Nh = Nr * (1 + 10 * abs(nu_r(6))) * nu_r(6);

tau_damp = [Xh Yh Zh Kh Mh Nh]';

% Strip theory: cross-flow drag integrals
tau_crossflow = crossFlowDrag(L,B_pont,T,nu_r);

% Payload expressed in NED
f_payload = Rzyx(eta(4),eta(5),eta(6))' * [ 0 0 mP*g ]';  % Payload force 
m_payload = Smtrx(r_bP) * f_payload;                        % Payload moment 
g_0 = [ f_payload
        m_payload ];

% Trim condition: G * eta_0 = g_0
eta_0 = [0; 0; inv(G(3:5,3:5)) * g_0(3:5); 0];
eta = eta - eta_0; % shifted equilibrium

% Kinematic transformation matrix
J = eulerang(eta(4),eta(5),eta(6));

% Time derivative of the state vector, numerical integration see ExOtter.m  
xdot = [ nu_c_dot + ...
         M \ ( tau + tau_damp + tau_crossflow - C * nu_r - G * eta )
         J * nu ];  
     
end
