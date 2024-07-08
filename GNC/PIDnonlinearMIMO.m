function tau = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h)
% tau = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h) 
% Nonlinear MIMO PID regulator for dynamic positioning (DP). For 3-DOF
% models (surge, sway, and yaw) the output is: tau = [tau1, tau2, tau6]'. 
% In the general case, tau = [tau1, tau2, 0, 0, 0, tau6]', where
%    .
%    z_int = eta - eta_d,     where eta_d = 1 / (T_f * s + 1) * eta_ref
%
%    tau = -R(psi)' * ( Kp * (eta - eta_d) + Ki * z_int ) - Kd * nu
%
%    Kp = M_diag * wn * wn,                  M_diag = diag(diag(M))
%    Kd = M_diag * 2 * zeta * wn
%    Ki = 1/10 * Kp * wn
%
% is based on Algorithm 15.2, MIMO nonlinear PID Pole-Placement Algorithm, 
% by Fossen (2021); see also Equation (15.82). 
%
% Persistent variables: 
% The setpoint eta_d and integral state z_int are persistent variables that 
% should be cleared by adding:
%
%    clear PIDnonlinearMIMO
%
% on the top in the script calling PIDnonlinearMIMO.m.
%
% Inputs:  
%   eta: generalized position vector, 3x1 (surge, sway, and yaw) or 6x1
%   nu:  generalized velocity vector, 3x1 (surge, sway, and yaw) or 6x1
%   eta_ref: vector [x_ref,y_ref,psi_ref] of setpoints in surge, sway, and yaw
%   M: system inertia matrix, 3x3 (surge, sway, and yaw) or 6x6
%   wn: closed-loop natural frequencies, scalar or diagonal matrix 3x3
%   zeta: closed-loop relative damping ratios, scalar or diagonal matrix 3x3
%   T_f: setpointlow-pass filter time constant (s)
%   h: sampling time (s)
%
% Outputs:  
%    tau: generalized control force, 3x1 (surge, sway, and yaw) or 6x1
%  
% Author:    Thor I. Fossen
% Date:      2 Sep 2023

persistent eta_d;  % LP-filtered commands
persistent z_int;  % integral states

% Initialization of desired state eta_d and integral state z_int 
if isempty(z_int)
    z_int = [0 0 0]';             
    eta_d = [0 0 0]';
end

eta_ref = eta_ref(:);  % Make sure eta_ref is a column vector

% Reduce 6-DOF model to 3-DOF model
DOF = 3;
if length(nu) == 6
    eta = [eta(1) eta(2) eta(6)]';
    nu  = [nu(1) nu(2) nu(6)]';
    M = M([1 2 6],[1 2 6]);
    DOF = 6;
end

% Rotation matrix in yaw
R = [ cos(eta(3)) -sin(eta(3)) 0
      sin(eta(3))  cos(eta(3)) 0
      0            0           1 ];

% MIMO pole placement (Algorithm 15.2)
M_diag = diag(diag(M));             % diagonal M matrix
Kp = M_diag .* wn .* wn;
Kd = M_diag .* 2 .* zeta .* wn;
Ki = 1/10 * Kp .* wn;

% 3-DOF control law for surge, sway and yaw
e = eta - eta_d;
e(3) = ssa( e(3) );
tau_PID = -R' * ( Kp * e + Ki * z_int ) - Kd * nu;

if DOF == 6
    tau = [tau_PID(1), tau_PID(2), 0, 0, 0, tau_PID(3)]'; 
else  % 3-DOF
    tau = tau_PID;
end

% integral state: z_int[k+1]
z_int = z_int + h * (eta - eta_d);

% low-pass filter: eta_d[k+1]
eta_d = eta_d + h * (eta_ref - eta_d)/ T_f;

end

