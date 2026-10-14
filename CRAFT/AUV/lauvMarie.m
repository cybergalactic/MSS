function xdot = lauvMarie(x,tau_ext,Vc,betaVc,w_c,w0_spectrum)
% lauvMarie  Unified 6-DOF seakeeping and maneuvering model for LAUV Marie.
% xdot = lauvMarie(x,tau_ext,Vc,betaVc,w_c,w0_spectrum) returns the
% derivative of a 12-state, 6-DOF model expressed in forward-starboard-down
% (FSD) body axes and North-East-Down (NED) navigation coordinates.
%
% The model combines seakeeping hydrodynamics and nonlinear maneuvering
% forces in one time-domain equation. Frequency-dependent added mass and
% potential damping from the MSS-Capytaine LAUV_marie data set are reduced
% to constant, spectrum-dependent matrices using the power-based method:
%
%     M = MRB + A_eq,                 D = B_eq + Bv
%
% A_eq and B_eq are weighted by a normalized wave spectrum with peak
% angular frequency w0_spectrum. The reduction preserves the spectrum-
% weighted kinetic-energy and dissipation properties of the original
% frequency-dependent model. This supplies a maneuvering model that remains
% consistent with the selected seakeeping condition without introducing a
% separate radiation-memory state model.
%
% The relative-velocity equations have the form:
%
%     eta_dot = J(eta) * nu
%
%     M * nu_r_dot + (CRB + CA) * nu_r + D(nu_r) * nu_r + g(eta)
%         = tau_ext + tau_crossflow + tau_liftdrag
%
% with the current-acceleration correction included when converting from
% nu_r_dot to nu_dot. Here nu_r is the body velocity relative to the ocean
% current. tau_ext is supplied by the caller and can contain actuator,
% control, and first-order wave-excitation forces and moments. Wave-force
% RAO evaluation, directional spreading, and depth attenuation are handled
% outside this function, for example by SIMlauvMarie.
%
% POWER-BASED AND VISCOUS DAMPING
%   The constant potential-flow matrices A_eq and B_eq are computed by
%   computeManeuveringModel. The additional diagonal viscous matrix Bv uses
%   target aperiodic time constants T_1236 for surge, sway, heave, and yaw,
%   together with damping-ratio increments delta_zeta_45 for roll and pitch.
%   The current nominal values are:
%
%       T_1236        = [50 5 5 5] seconds
%       delta_zeta_45 = [0.2 0.2]
%
%   These quantities are engineering model parameters, not immutable vehicle
%   properties. They can be identified or calibrated using captive tests,
%   free-decay tests, towing-tank measurements, or full-scale input-output
%   data. Calibration should cover the speeds, depths, and sea states for
%   which the model will be used.
%
% NONLINEAR MANEUVERING TERMS
%   The constant power-based model is augmented by:
%     - exponential blending of linear surge damping at increasing speed;
%     - quadratic roll, pitch, and yaw damping multipliers;
%     - strip-wise cross-flow drag using a circular-cylinder section model;
%     - angle-of-attack-dependent lift and drag; and
%     - nonlinear rigid-body and added-mass Coriolis-centripetal terms.
%
%   The nonlinear damping multipliers, drag coefficients, and lift/drag
%   parameters are also calibratable. The supplied values are baseline
%   assumptions for the idealized LAUV geometry and should be validated
%   against experiments before quantitative prediction or controller tuning.
%
% MODELING ASSUMPTIONS
%   The hydrodynamic data represent an idealized, bare LAUV Marie hull and do
%   not include fins, antennas, sonar heads, or the propeller. Selected
%   added-mass Coriolis couplings are suppressed because retaining them
%   without the missing stabilizing appendage derivatives would produce an
%   unbalanced Munk-moment model. Neutral buoyancy is assumed, with restoring
%   moments produced by the specified CG-CB separation. Actuator forces,
%   actuator dynamics, and first-order wave loads must be included through
%   tau_ext. Consequently, this is a unified reduced-order seakeeping and
%   maneuvering model, but it is not a validated digital twin of the vehicle.
%
% INPUTS
%   x             12x1 state vector
%                 [u v w p q r north east down phi theta psi]'. Linear and
%                 angular velocities are in m/s and rad/s; NED positions and
%                 Euler angles are in m and rad.
%   tau_ext       6x1 external generalized force/moment vector
%                 [X Y Z K M N]' in N and N m. It normally contains the sum
%                 of control/propulsion loads and first-order wave loads.
%   Vc            Horizontal ocean-current speed (m/s).
%   betaVc        Horizontal current direction in NED coordinates (rad).
%   w_c           Vertical ocean-current speed, positive down (m/s).
%   w0_spectrum   Wave-spectrum peak angular frequency (rad/s) used for the
%                 power-based reduction of added mass and potential damping.
%
% OUTPUT
%   xdot          12x1 state derivative [nu_dot; eta_dot]
%
% PERSISTENT INITIALIZATION
%   The vessel data and power-based matrices are initialized on the first
%   call and then cached. If w0_spectrum or any damping calibration parameter
%   is changed between simulations, execute "clear lauvMarie" before the next
%   call so that the equivalent model is recomputed.
%
% REFERENCES
%   Fossen, T. I. (2025). Maneuvering Coefficient Estimation from Frequency-
%   Dependent Added Mass and Damping: A Power-Based Approach. Ocean
%   Engineering, 341, 122494.
%
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-10-14
% Revisions:

persistent vessel

% If the 'vessel' strcuture is missing, the following quantities are computed once
if isempty(vessel)
    matFile = resolveHydroVesselFile('LAUV_marie.mat');
    auvData = load(matFile, 'vessel');
    vessel = auvData.vessel;

    % Power-based equivalent matrices and viscous damping
    vessel.powerBased.T_1236 = [50 5 5 5];
    vessel.powerBased.delta_zeta_45 = [0.2 0.2];

    vessel = computeManeuveringModel(vessel, w0_spectrum, ...
        vessel.powerBased.T_1236, vessel.powerBased.delta_zeta_45, 0);

    % Nonlinear damping multipliers (quadratic damping) in roll, pitch, and yaw
    vessel.kappa_4 = 5;
    vessel.kappa_5 = 5;
    vessel.kappa_6 = 0;
end

% MSS global constants
constant = mssConstants();

% Decompose the state into body-fixed velocity and NED position/attitude
nu = x(1:6);       % Generalized velocity [u v w p q r]'
eta = x(7:12);     % Position and Euler angles [N E D phi theta psi]'

% Transform the irrotational North-East current into body-fixed coordinates.
% Its body-fixed derivative is caused solely by the vessel angular velocity.
v_c = [ Vc * cos(betaVc - eta(6))
    Vc * sin(betaVc - eta(6))
    w_c ];
nu_c = [v_c; zeros(3,1)];
nu_c_dot = [-Smtrx(nu(4:6)) * v_c; zeros(3,1)];
nu_r = nu - nu_c;  % Vessel velocity relative to the current

% Coriolis-centripetal matrices. The rigid-body representation is independent
% of the relative linear velocity [u_r, v_r, w_r]'.
[~, CRB] = rbody(vessel.main.m, ...
    vessel.main.k44, vessel.main.k55, vessel.main.k66, nu(4:6), [0 0 0]');
CA = m2c(vessel.MA, nu_r);

% Augment the linear damping matrix with nonlinear diagonal damping. The
% lift-drag model below supplies the quadratic surge resistance, while the
% exponential term blends out the linear surge damping at higher speed.
k_u = 3;  % Exponential blending coefficient (s/m)
D_nonlinear = vessel.D;
D_nonlinear(1,1) = vessel.D(1,1) * exp(-k_u * abs(nu_r(1)));
D_nonlinear(4,4) = vessel.D(4,4) * (1 + vessel.kappa_4 * abs(nu(4)));
D_nonlinear(5,5) = vessel.D(5,5) * (1 + vessel.kappa_5 * abs(nu(5)));
D_nonlinear(6,6) = vessel.D(6,6) * (1 + vessel.kappa_6 * abs(nu(6)));

% The MSS-Capytaine LAUV model omits stabilizing fins, appendages, and
% validated maneuvering derivatives. The selected added-mass Coriolis
% couplings are therefore omitted to avoid retaining an unbalanced
% Munk-moment model. This is a reduced-order modeling assumption; the
% corresponding physical terms are not generally zero for a bare body.
CA(5,3) = 0; CA(3,5) = 0;  % Heave-pitch coupling
CA(5,1) = 0; CA(1,5) = 0;  % Surge-pitch Munk coupling
CA(6,1) = 0; CA(1,6) = 0;  % Yaw-related Munk couplings
CA(6,2) = 0; CA(2,6) = 0;

% Using cross-flow drag for cylinders and adding lift-drag forces
tau_crossflow = crossFlowDrag(...
    vessel.main.Lpp, vessel.main.B, vessel.main.T, nu_r, 'cylinder');

L_auv = vessel.main.Lpp;
D_auv = vessel.main.T;
R_auv = D_auv / 2;
S = 0.7 * L_auv * D_auv; % Planform area S = 70% of rectangle L_auv * D_auv
Cd = 0.42; % From Allen et al. (2000)
CD_0 = Cd * pi * R_auv^2 / S; % Parasitic drag coefficient CD_0 (alpha = 0)
alpha = atan2( nu_r(3), nu_r(1) ); % Angle of attack (rad)
U_r = norm(nu_r(1:3)); % Relative speed (m/s)
tau_liftdrag = forceLiftDrag(D_auv, S, CD_0, alpha, U_r);

% Restoring forces and moments
W = vessel.main.m * constant.g;  B = W;
g_eta = gvect(W, B, eta(5), eta(4), vessel.main.CG, vessel.main.CB);

% 6-DOF kinetics in body-fixed FSD coordinates
nudot = nu_c_dot + vessel.Minv * (tau_ext + tau_crossflow + tau_liftdrag ...
    - (CRB + CA + D_nonlinear) * nu_r - g_eta);

% Transform body-fixed velocity to NED position and Euler-angle rates
J = eulerang(eta(4),eta(5),eta(6));
etadot = J * nu;

% Assemble the state derivative for numerical integration
xdot = [nudot; etadot];

end
