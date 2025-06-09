function tau_crossflow = crossFlowDrag(L,B,T,nu_r,drag_model)
% tau_crossflow = crossFlowDrag(L,B,T,nu_r) computes the cross-flow drag 
% integrals for a marine craft using strip theory. Application:
%
%  M d/dt nu_r + C(nu_r)*nu_r + D*nu_r + g(eta) = tau + tau_crossflow
%
% Inputs: L:  length
%         B:  beam
%         T:  draft 
%         nu_r = [u-u_c, v-v_c, w-w_c, p, q, r]': relative velocity vector
%
% Output: tau_crossflow = [0 Yh Zh 0 Mh Nh]: 6-DOF cross-flow drag vector
%
% Author:     Thor I. Fossen 
% Date:       25 Apr 2021, Horizontal-plane drag of ships
% Revisions:  30 Jan 2021, Extended to include heave and pitch for AUVs
%             09 Jun 2025, Make different drag models selectable (M. Seidl)

if (nargin == 4), drag_model = 'Hoerner'; end % Hoerner is default drag model

rho = 1025;             % density of water

dx = L / 20;            % divide marine craft into 20 strips

switch drag_model
    case 'Hoerner'
        Cd_2D = Hoerner(B,T);   % 2-D drag coefficient based on Hoerner's curve
    case 'cylinder'
        Cd_2D = cylinderDrag(L,B,nu_r); % 2D drag coefficient based on cylinder data from DNV-RP-C205
    otherwise
        error('Unsupported drag model %s.',drag_model)
end

Yh = 0; Zh = 0; Mh = 0; Nh = 0;
for xL = -L/2:dx:L/2
    v_r = nu_r(2);          % relative sway velocity
    w_r = nu_r(3);
    q = nu_r(5);            % pitch rate
    r = nu_r(6);            % yaw rate
    U_h = abs(v_r + xL * r) * (v_r + xL * r);
    U_v = abs(w_r + xL * q) * (w_r + xL * q);   
    Yh = Yh - 0.5 * rho * T * Cd_2D * U_h * dx;       % sway force
    Zh = Zh - 0.5 * rho * T * Cd_2D * U_v * dx;       % heave force  
    Mh = Mh - 0.5 * rho * T * Cd_2D * xL * U_v * dx;  % pitch moment    
    Nh = Nh - 0.5 * rho * T * Cd_2D * xL * U_h * dx;  % yaw moment
end

tau_crossflow = [0 Yh Zh 0 Mh Nh]';

end

