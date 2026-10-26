function [xdot,U,M,B_prop,n_min,n_max] = otter(x,n,mP,r_bP,V_c,beta_c)
% Unified 6-DOF maneuvering and seakeeping dynamics for the Otter USV.
% [xdot,U,M,B_prop,n_min,n_max] = otter(x,n,mP,r_bP,V_c,beta_c) returns the
% time derivative xdot of the state vector, the speed U in m/s, the 6x6
% mass matrix M, the 2x2 propeller input matrix B_prop, and the physical
% propeller shaft-speed limits n_min and n_max.
%
%  x = [ u v w p q r x y z phi theta psi ]'
%
% The model describes the 2.0 m Maritime Robotics Otter USV; see
% www.maritimerobotics.com. It combines rigid-body and payload dynamics,
% hydrostatic restoring, propeller thrust, ocean currents, and nonlinear
% damping with added mass and radiation damping computed by MSS-Capytaine.
% Nonlinear rigid-body and added-mass Coriolis/centripetal terms are included,
% with the selected pitch/yaw couplings omitted as documented below. SIMotter
% completes the unified simulation with first-order Capytaine force-RAO loads.
%
% Hydrodynamic model:
%   otter_80kg.mat contains frequency-dependent added mass A(omega),
%   potential damping B(omega), and wave RAO data for a 55 kg hull with a
%   nominal 25 kg payload. computeManeuveringModel uses normalized
%   wave-spectrum weighting to obtain constant power-based equivalent
%   matrices A_eq and B_eq. A fixed reference spectrum is used; changing the
%   SIMotter wave inputs does not recompute these cached matrices.
%   Additional viscous damping Bv uses surge, sway and yaw decay time
%   constants and heave, roll and pitch damping-ratio increments.
%
%   The matrices are referenced to the fixed CG stored in the MAT-file.
%   They are transformed once to the body-frame origin CO and cached:
%       MA = H' * A_eq * H
%       D  = H' * (B_eq + Bv) * H
%       G  = H' * C * H
%   D remains constant; D_nonlinear adds velocity-dependent damping at
%   each evaluation. Cross-flow drag is summed separately over both
%   pontoons using local velocities and forces/moments transformed to CO.
%   The strip model treats the pontoons independently and neglects their
%   viscous interaction; lateral reference offsets are [0, +/-y_pont, 0].
%
%   Rigid-body inertia is recomputed for the current payload mass and
%   position. CRB is evaluated at the current angular velocity.
%   Payload gravity shifts the restoring equilibrium; the linear
%   stiffness G, hydrodynamic data, and viscous damping calibration remain
%   based on the nominal 80 kg load condition. Damping increments and nonlinear
%   factors are tuning parameters that require experimental calibration.
%
% Wave loads:
%   SIMotter applies first-order excitation from vessel.forceRAO to either
%   dynamics model. It transforms forces from the data reference CG to CO
%   and adds M\tau_wave to the acceleration at each RK4 stage. The present
%   function retains its original calling interface and includes equivalent
%   radiation coefficients; external wave excitation is applied by SIMotter.
%
% otter_legacy.m retains the previous first-principles implementation.
% The state vector and calling interface are retained for existing callers.
%
% The state vector is defined as:
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
% nominal draft is T = 13.4 cm for mP = 0 and T = 19.5 cm for mP = 25 kg.
% The unloaded condition defines z = 0. Payload gravity is included through
% g_0 and shifts the restoring equilibrium. Kinematics use the actual
% attitude, while restoring forces use displacement from this equilibrium.
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
%   SIMotter.m : Heading/path-following simulation with currents and waves.
%   exOtter.m  : Legacy-model example demonstrating the 5-state EKF for SOG,
%                COG, and course-rate estimation with a course autopilot.
%
% References:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%     Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
%   Fossen, T. I. (2025). Maneuvering Coefficient Estimation from Frequency-
%     Dependent Added Mass and Damping: A Power-Based Approach.
%     Ocean Engineering, 341, 122494.
%
% Author:    Thor I. Fossen
% Date:      2026-10-26

% The constant hydrodynamic matrices MA and D are computed only once
persistent MA D G;

% ------------------------------------------------------------------------------
% Check of inputs and state dimensions
% ------------------------------------------------------------------------------
if nargin == 0
    x = zeros(12,1); n = zeros(2,1); mP=25; r_bP = zeros(3,1); V_c=0; beta_c=0;
end

if nargin < 5
    V_c = 0;
end

if nargin < 6
    beta_c = 0;
end

if (length(x) ~= 12),error('x vector must have dimension 12!'); end
if (length(n) ~= 2),error('n vector must have dimension 2!'); end

% ------------------------------------------------------------------------------
% Main data
% ------------------------------------------------------------------------------
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
zeta_heave = 0;       % Additional heave damping ratio, potential damping B_eq(3,3) is small
zeta_roll = 0.2;      % Additional roll damping ratio, potential damping B_eq(4,4) is small
zeta_pitch = 0.1;     % Additional pitch damping ratio, potential damping B_eq(5,5) is small
kappa_4 = 5;          % Quadratic roll damping factor
kappa_5 = 0;          % Quadratic pitch damping factor
kappa_6 = 5;          % Quadratic yaw damping factor
Umax = 6 * 0.5144;    % 6 knots maximum forward speed, used to compute Xuu at max thrust

% Compute hydrodynamic data only once to avoid that RK4 repeats the computation
% in the loop at each time step
if isempty(MA)

    % Load MSS-Capytaine hydrodynamic data for 55 kg Otter with 25 kg payload
    load otter_80kg vessel;

    % Compute power-based equivalent added mass/damping and viscous damping
    vessel.powerBased.T_126 = [T_surge T_sway T_yaw];
    vessel.powerBased.delta_zeta_345 =[zeta_heave zeta_roll zeta_pitch];
    omega_p = 3.0;  % Fixed reference spectrum for equivalent matrices
    vessel = computeManeuveringModel(vessel, omega_p);

    % The Capytaine hydrodynamic data are at the fixed CG stored in the MAT-file.
    % Transform once to CO; the payload-dependent CG below applies only
    % to the rigid-body matrices, not to this hydrodynamic reference point.
    H_hydro = Hmtrx(vessel.main.CG(:));
    MA = H_hydro' * vessel.powerBased.A_eq * H_hydro;
    D = H_hydro' * (vessel.powerBased.Bv + vessel.powerBased.B_eq) * H_hydro;

    % Capytaine hydrostatics include the heave-pitch coupling associated
    % with vessel.main.CF(1) = 0.017 m for this 80 kg dataset.
    G = H_hydro' * vessel.C(:,:,1) * H_hydro;
end

% Data for one pontoon
B_pont  = 0.25;       % Beam of one pontoon (m)
y_pont  = 0.395;      % Distance from centerline to waterline area center (m)
Cb_pont = 0.4;        % Block coefficient, computed from m = 55 kg

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

% Transform rigid-body matrices from the current CG to CO
H = Hmtrx(r_bG);
MRB = H' * MRB_CG * H;
CRB = H' * CRB_CG * H;

% Coriolis and centripetal matrix
CA = m2c(MA, nu_r);

% The Otter model does not retain a complete set of nonlinear maneuvering
% derivatives for the hull, thrusters, and appendages. Selected added-mass
% Coriolis couplings in pitch and yaw are therefore omitted to avoid retaining
% unbalanced Munk-moment contributions without the complementary terms.
% For a full added-mass matrix, these entries include additional couplings.
% This is a modeling approximation; the physical terms need not vanish.
% Keep the full MA and cancel entries in pairs to preserve CA skew symmetry.
CA(5,3) = 0; CA(3,5) = 0;  % Heave-pitch coupling
CA(5,1) = 0; CA(1,5) = 0;  % Surge-pitch coupling
CA(6,1) = 0; CA(1,6) = 0;  % Surge-yaw coupling
CA(6,2) = 0; CA(2,6) = 0;  % Sway-yaw coupling

% System mass and Coriolis-centripetal matrices
M = MRB + MA;
C = CRB + CA;

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

% Nonlinear damping matrix
D_nonlinear = D;

% Blending linear and quadratic damping
Xuu = -24.4 * g / Umax^2;  % Specified using the maximum speed
D_nonlinear(1,1) = D(1,1) * exp(-3 * abs(nu_r(1))) - Xuu * abs(nu_r(1));

% Add quadratic roll, pitch and yaw damping
D_nonlinear(4,4) = D(4,4) * (1 + kappa_4 * abs(nu_r(4)));
D_nonlinear(5,5) = D(5,5) * (1 + kappa_5 * abs(nu_r(5)));
D_nonlinear(6,6) = D(6,6) * (1 + kappa_6 * abs(nu_r(6)));

% Strip theory: sum independent pontoon drag contributions at CO
tau_crossflow = zeros(6,1);
for y = [-y_pont, y_pont]
    H_pont = Hmtrx([0 y 0]);
    nu_pont = H_pont * nu_r;
    tau_pont = crossFlowDrag(L,B_pont,T,nu_pont);
    tau_crossflow = tau_crossflow + H_pont' * tau_pont;
end

% Payload expressed in NED
f_payload = Rzyx(eta(4),eta(5),eta(6))' * [ 0 0 mP*g ]';  % Payload force
m_payload = Smtrx(r_bP) * f_payload;                      % Payload moment
g_0 = [ f_payload
        m_payload ];

% Trim condition: G * eta_0 = g_0
eta_0 = [0; 0; G(3:5,3:5) \ g_0(3:5); 0];

% Kinematic transformation matrix
J = eulerang(eta(4),eta(5),eta(6));

% Time derivative of the state vector
xdot = [ nu_c_dot + ...
         M \ (tau + tau_crossflow - D_nonlinear * nu_r - C * nu_r - G * (eta-eta_0))
         J * nu ];

end
