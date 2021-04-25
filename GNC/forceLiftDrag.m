function tau_liftdrag = forceLiftDrag(b,S,CD_0,alpha,U_r)
% tau_liftdrag = forceLiftDrag(b,S,CD_0,alpha,Ur) computes the hydrodynamic
% lift and drag forces of a submerged "wing profile" for varying angle of
% attack (Beard and McLain 2012). Application:
%
%  M d/dt nu_r + C(nu_r)*nu_r + D*nu_r + g(eta) = tau + tau_liftdrag
%
% Output:
%  tau_liftdrag:  6x1 generalized force vector
%
% Inputs:
%  b:       wing span (m)
%  S:       wing area (m^2)
%  CD_0:    parasitic drag (alpha = 0), typically 0.1-0.2 for a streamlined body
%  alpha:   angle of attack, scalar or vector (rad)
%  U_r:     relative speed (m/s)
%
% Example:
%
% Cylinder-shaped AUV with length L = 1.8, diameter D = 0.2 and CD_0 = 0.1:
%    tau_liftdrag = forceLiftDrag(0.2, 1.8*0.2, 0.1, alpha, U_r)
% 
% Author:    Thor I. Fossen
% Date:      25 April 2021 

rho = 1026;

[CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,0);

F_drag = 1/2 * rho * U_r^2 * S * CD;    % drag force
F_lift = 1/2 * rho * U_r^2 * S * CL;    % lift force

% transform from FLOW axes to BODY axes using angle of attack
tau_liftdrag = [...
    cos(alpha) * (-F_drag) - sin(alpha) * (-F_lift)
    0
    sin(alpha) * (-F_drag) + cos(alpha) * (-F_lift)
    0
    0
    0 ];