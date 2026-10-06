function xdot = hydroVessel(x, tau_ext, vessel, Vc, betaVc)
% xdot = hydroVessel(x,tau_ext,vessel,Vc,betaVc) returns the time derivative
% of a 12-state, 6-DOF vessel model. Body-fixed quantities use the
% front-starboard-down (FSD) convention, and all matrices and generalized
% forces are computed about the center of gravity (CG).
%
% The equations of motion are
%
%   eta_dot = J(eta) * nu
%   nu_dot  = nu_c_dot + Minv * (tau_ext + tau_crossflow
%              - (CRB + CA + D_nonlinear) * nu_r - G * eta)
%
% where nu_r = nu - nu_c is the velocity relative to the ocean current.
%
% Inputs:
%   x:        12x1 state vector
%             x = [u v w p q r north east down phi theta psi]'
%             Linear and angular velocities are in m/s and rad/s; positions
%             and Euler angles are in m and rad, respectively.
%   tau_ext:  6x1 external generalized force vector [N; N; N; Nm; Nm; Nm].
%             In SIMhydroVessel, tau_ext = tau + tau_wave1.
%   vessel:   Structure containing the constant maneuvering model:
%               vessel.MA, vessel.Minv  - added-mass and inverse mass matrices
%               vessel.D, vessel.G      - damping and restoring matrices
%               vessel.kappa_4,5,6      - nonlinear damping coefficients
%               vessel.main             - mass, radii of gyration, density,
%                                         dimensions, and hull data
%   Vc:       Horizontal ocean-current speed (m/s)
%   betaVc:   Ocean-current direction in the North-East frame (rad)
%
% Output:
%   xdot:     12x1 state derivative [nu_dot; eta_dot]
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2026-10-02
% Revisions:

% Decompose the state into body-fixed velocity and NED position/attitude
nu = x(1:6);       % Generalized velocity [u v w p q r]'
eta = x(7:12);     % Position and Euler angles [N E D phi theta psi]'

% Transform the irrotational North-East current into body-fixed coordinates.
% Its body-fixed derivative is caused solely by the vessel angular velocity.
v_c = [ Vc * cos(betaVc - eta(6))
        Vc * sin(betaVc - eta(6))
                                0 ];
nu_c = [v_c; zeros(3,1)];
nu_c_dot = [-Smtrx(nu(4:6)) * v_c; zeros(3,1)];
nu_r = nu - nu_c;  % Vessel velocity relative to the current

% Coriolis-centripetal matrices. The rigid-body representation is independent
% of the relative linear velocity [u_r, v_r, w_r]'.
[~, CRB] = rbody(vessel.main.m, ...
    vessel.main.k44, vessel.main.k55, vessel.main.k66, nu(4:6), [0 0 0]');
CA = m2c(vessel.MA, nu_r);

% ITTC quadratic surge drag with exponential blending of the linear damping.
if strcmpi(vessel.main.name,'semisub')

    % ITTC quadratic surge drag is not used for the low-speed semisubmersible
    Xuu = 0;
    k_u = 0;  % Retain the full linear surge damping
else

    % Ships and USVs
    if isfield(vessel.main,'C_B') && vessel.main.C_B > 0
        C_B = vessel.main.C_B;
    elseif isfield(vessel.main,'nabla') && vessel.main.nabla > 0
        C_B = vessel.main.nabla / ...
            (vessel.main.Lpp * vessel.main.B * vessel.main.T);
        C_B = min(C_B,1);
    else
        C_B = [];  % Use the default value in XuuITTC
    end

    Xuu = XuuITTC(nu_r(1),vessel.main.rho, ...
        vessel.main.Lpp,vessel.main.B,vessel.main.T,C_B);
    k_u = 3;  % Exponential blending coefficient (s/m)
end

% Augment the linear damping matrix with nonlinear diagonal damping. Surge
% blends from linear damping at low speed to the ITTC quadratic resistance;
% roll, pitch, and yaw use user-defined quadratic damping multipliers.
D_nonlinear = vessel.D;
D_nonlinear(1,1) = vessel.D(1,1) * exp(-k_u * abs(nu_r(1))) - Xuu * abs(nu_r(1));
D_nonlinear(4,4) = vessel.D(4,4) * (1 + vessel.kappa_4 * abs(nu(4)));
D_nonlinear(5,5) = vessel.D(5,5) * (1 + vessel.kappa_5 * abs(nu(5)));
D_nonlinear(6,6) = vessel.D(6,6) * (1 + vessel.kappa_6 * abs(nu(6)));

% Add cross-flow drag and lift-drag forces
tau_liftdrag = zeros(6,1); % Zero for displacement vessels
if strcmpi(vessel.main.name,'semisub')
    % Using the dimensions of each pontoon from WAMIT geometry semisub.gdf
    L_pontoon = 115; B_pontoon = 18; T_pontoon = 10;
    tau_crossflow = 2 * crossFlowDrag(...
        L_pontoon, B_pontoon, T_pontoon, nu_r, 'Hoerner');
elseif strcmpi(vessel.main.name,'LAUV_marie')

    % The idealized LAUV model omits stabilizing fins, appendages, and
    % validated maneuvering derivatives. The selected added-mass Coriolis
    % couplings are therefore omitted to avoid retaining an unbalanced
    % Munk-moment model. This is a reduced-order modeling assumption; the
    % corresponding physical terms are not generally zero for a bare body.
    CA(5,3) = 0; CA(3,5) = 0;  % Heave-pitch coupling
    CA(5,1) = 0; CA(1,5) = 0;  % Surge-pitch Munk coupling
    CA(6,1) = 0; CA(1,6) = 0;  % Yaw-related Munk couplings
    CA(6,2) = 0; CA(2,6) = 0;

    % Using cross-flow drag for cylinders and adding lift-drag forces
    D_nonlinear(1,1) = vessel.D(1,1) * exp(-k_u * abs(nu_r(1)));
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
else
    % Using Hoerner's curve for monohull cross-flow drag
    tau_crossflow = crossFlowDrag(...
        vessel.main.Lpp, vessel.main.B, vessel.main.T, nu_r, 'Hoerner');
end

% 6-DOF kinetics in body-fixed FSD coordinates
nudot = nu_c_dot + vessel.Minv * (tau_ext + tau_crossflow + tau_liftdrag ...
    - (CRB + CA + D_nonlinear) * nu_r - vessel.G * eta);

% Transform body-fixed velocity to NED position and Euler-angle rates
J = eulerang(eta(4),eta(5),eta(6));
etadot = J * nu;

% Assemble the state derivative for numerical integration
xdot = [nudot; etadot];

end
