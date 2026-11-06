function delta = integralSMCheading(...
    psi,r,psi_d,r_d,a_d,K_d,K_sigma,lambda,phi_b,K_nomoto,T_nomoto,h,delta_max)
% PID and integral sliding mode controller (SMC) for heading control. The 
% yaw dynamics is modelled by 
%
%  psi_dot = r
%  r_dot + (1 / T_yaw) * r = (K_nomoto / T_nomoto) * delta
% 
% where delta is the rudder angle. The input gain, K_nomoto, is found from 
% the steady-state condition, r_max = K_nomoto * delta_max, while T_nomoto 
% is the time constant in yaw observered during steady turning.
%
% The heading autopilot (Equation 16.479 in Fossen 2027) sliding surface
% and control law are
%
%   sigma = r-r_d + 2*lambda*ssa(psi-psi_d) + lambda^2 * integral(ssa(psi-psi_d))
%
%   delta = (T_nomoto * r_r_dot + r_r - K_d * sigma 
%       - K_sigma*(sigma/phi_b)) / K_nomoto
% 
% where lambda > 0, K_d > 0 (PID control), and K_sigma > 0 (SMC). 
% The integral state psi_int is a persistent variable that should be cleared 
% by adding:
%
%    clear integralSMCheading
%
% on the top in the script calling integralSMCheading.m.
%
% Inputs:  
%   psi: yaw angle (rad)
%   r: yaw rate(rad/s)
%   psi_d: desired yaw angle (rad)
%   r_d: desired yaw rate(rad/s)
%   a_d: desired yaw acceleration (rad/s^2)
%   K_d: PID control gain, multiplied with the sliding surface sigma
%   K_sigma: SMC gain
%   lambda: constant defining the dynamics on the slideing surface
%   phi_b: boundary layer thickness
%   K_nomto: Nomoto gain constant (1/s)
%   T_nomto: Nomoto time constant (s)
%   h: sampling time (s)
%   delta_max: (optional) rudder angle limit (rad). When given, the output is
%              limited to +-delta_max and the integral state is not updated
%              while the limit is active (conditional integration, no integral
%              windup). Without it, the output is not limited.
%
% Outputs:  
%    delta: rudder angle
%  
% Author:    Thor I. Fossen
% Date:      2024-02-09
% Revisions:
%   2026-11-06: Optional rudder limit delta_max with conditional integration;
%               the default is unchanged.

persistent psi_int;               % integral state

% Initialization of desired state psi_d and integral state z_int 
if isempty(psi_int)
    psi_int = 0;             
end

% PID and integral SMC (Equation 16.479 in Fossen 2027)
e_psi = ssa( psi - psi_d );
e_r = r - r_d;

r_r_dot = a_d - 2 * lambda * e_r - lambda^2 * e_psi;
r_r     = r_d - 2 * lambda * e_psi - lambda^2 * psi_int;
sigma   = r - r_r;

if abs(sigma / phi_b) > 1
    delta = ( T_nomoto * r_r_dot + r_r - K_d * sigma...
        - K_sigma * sign(sigma) ) / K_nomoto;
else
    delta = ( T_nomoto * r_r_dot + r_r - K_d * sigma...
        - K_sigma * (sigma / phi_b) ) / K_nomoto;
end

% Optional output limit with conditional integration
limited = false;
if nargin > 12 && ~isempty(delta_max)
    if abs(delta) > delta_max
        limited = true;
        delta = sign(delta) * delta_max;
    end
end

% Propagation of persistent integral state: psi_int[k+1]
if ~limited
    psi_int = psi_int + h * e_psi;
end

end
