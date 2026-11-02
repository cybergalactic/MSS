function [nu_c, nu_c_dot] = currentVelocity(Vc, betaVc, w_c, Theta, nu2, currentModel)
% [nu_c, nu_c_dot] = currentVelocity(Vc, betaVc, w_c, Theta, nu2, currentModel)
% computes the generalized ocean current velocity nu_c = [u_c v_c w_c 0 0 0]'
% expressed in BODY and its time derivative in BODY for an irrotational
% current that is constant in NED (Fossen 2011, Section 8.3):
%
%   v_c^n = [ Vc * cos(betaVc), Vc * sin(betaVc), w_c ]'
%   v_c^b = R(Theta)' * v_c^n
%   d/dt v_c^b = -S(nu2) * v_c^b
%
% The time derivative follows from d/dt R = R * S(nu2) (Fossen 2011, Eq.
% 2.33 and Eqs. 8.139-8.141). It is used in the equations of motion
% expressed in the relative velocity nu_r = nu - nu_c,
%
%   nu_dot = nu_c_dot + M \ (tau - C(nu_r) * nu_r - D(nu_r) * nu_r - g(eta))
%
% Inputs:
%   Vc:           horizontal current speed (m/s)
%   betaVc:       current direction in NED (rad)
%   w_c:          vertical current speed, positive down (m/s)
%   Theta:        Euler angles [phi theta psi]' (rad), or unit quaternion
%                 q = [eta eps1 eps2 eps3]'
%   nu2:          angular velocity [p q r]' (rad/s)
%   currentModel: (optional) 'full_attitude' (default),
%                 'heading_only_rotation_rate' or 'heading_only_yaw_rate'
%
% Current models:
%   'full_attitude':              v_c^b = R(Theta)' * v_c^n and
%                                 d/dt v_c^b = -S(nu2) * v_c^b.
%   'heading_only_rotation_rate': v_c^b = Rz(psi)' * v_c^n and
%                                 d/dt v_c^b = -S(nu2) * v_c^b.
%   'heading_only_yaw_rate':      v_c^b = Rz(psi)' * v_c^n and
%                                 d/dt v_c^b = [r * v_c, -r * u_c, 0]'.
%   The two heading-only models neglect the roll and pitch angles (and the
%   roll and pitch rates in the last model), as in the 3-DOF approximation
%   R(Theta) = Rz(psi) for small phi and theta (Fossen 2011, Eq. 2.47).
%   They are exact for phi = theta = 0 and p = q = 0. The first is the
%   current model used by otter.m, osv.m and hydroVessel.m before
%   2026-10-23, the second that of remus100.m and npsauv.m. All five
%   models now use 'full_attitude' by default and accept the setting as
%   their last input.
%
% Outputs:
%   nu_c:         6x1 generalized current velocity in BODY
%   nu_c_dot:     6x1 time derivative of nu_c in BODY
%
% Example:
%   [nu_c, nu_c_dot] = currentVelocity(0.5, deg2rad(30), 0, x(10:12), x(4:6))
%
% Reference:
%   T. I. Fossen (2011). Handbook of Marine Craft Hydrodynamics and Motion
%       Control. Wiley. Eqs. 2.33 (p. 25), 2.47 (p. 27) and 8.137-8.159
%       (pp. 222-225).
%
% Author:    Enio Krizman
% Date:      2026-10-23

if nargin < 6 || isempty(currentModel)
    currentModel = 'full_attitude';
end

switch currentModel
    case 'full_attitude'
        if length(Theta) == 3
            R = Rzyx(Theta(1), Theta(2), Theta(3));
        else
            R = Rquat(Theta);
        end
        v_c = R' * [ Vc * cos(betaVc); Vc * sin(betaVc); w_c ];
        nu_c = [ v_c; zeros(3,1) ];
        nu_c_dot = [ -Smtrx(nu2) * v_c; zeros(3,1) ];

    case {'heading_only_rotation_rate', 'heading_only_yaw_rate'}
        if length(Theta) == 3
            psi = Theta(3);
        else
            psi = atan2(2*(Theta(1)*Theta(4) + Theta(2)*Theta(3)), ...
                1 - 2*(Theta(3)^2 + Theta(4)^2));
        end
        u_c = Vc * cos( betaVc - psi );
        v_c = Vc * sin( betaVc - psi );
        nu_c = [ u_c v_c w_c 0 0 0 ]';
        if strcmp(currentModel, 'heading_only_rotation_rate')
            nu_c_dot = [ -Smtrx(nu2) * nu_c(1:3); zeros(3,1) ];
        else
            nu_c_dot = [ nu2(3)*v_c -nu2(3)*u_c 0 0 0 0 ]';
        end

    otherwise
        error(['currentModel must be ''full_attitude'', ', ...
            '''heading_only_rotation_rate'' or ''heading_only_yaw_rate''']);
end

end
