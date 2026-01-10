function [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bG)
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bG) computes the 6x6 rigid-body
% mass and Coriolis–centripetal matrices expressed at the body-fixed origin CO.
% If r_bG = [0 0 0]', the center of gravity (CG) coincides with CO and the
% matrices are computed about CO = CG.
%
% In general nu = [u, v, w, p, q, r]' but for a rigid body there exists a
% parametrization CRB = CRB(nu2) which is independent of the linear velocity
% nu1 = [u, v, w]'. This is particularly useful for relative equations of
% motion with irrotational current velocity nu_c = [u_c, v_c, w_c, 0, 0, 0]':
%
%   MRB*nudot + CRB(nu2)*nu = MRB*nudot_r + CRB(nu2_r)*nu_r = tau
%
% where nu_r = nu - nu_c and nu2_r is the angular part of nu_r.
%
% Outputs:
%  MRB            - Rigid-body mass matrix about CO
%  CRB = CRB(nu2) - Coriolis–centripetal matrix about CO (independent of nu1)
%
% Inputs:
%  m              - Mass
%  R44,R55,R66    - Radii of gyration about the CG
%  nu2 = [p,q,r]' - Angular velocity vector
%  r_bG           - Lever arm from CO to CG, expressed in {b}
%
% Author:    Thor I. Fossen
% Date:      2019-04-19
% Revisions:
%  2021-04-21 - Improved the documentation

I3 = eye(3);
O3 = zeros(3,3);

% MRB and CRB expressed about the CG
I_G   = m * diag([R44^2, R55^2, R66^2]);
MRB_CG = [ m*I3   O3
           O3     I_G ];

CRB_CG = [ m*Smtrx(nu2)     O3
           O3              -Smtrx(I_G*nu2) ];

% Transform MRB and CRB from CG to CO using the screw transformation
H   = Hmtrx(r_bG);
MRB = H' * MRB_CG * H;
CRB = H' * CRB_CG * H;