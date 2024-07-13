function [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bp)
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bp) computes the 6x6 rigid-body 
% mass and Coriolis-centripetal matrices for an arbitrariliy point P with
% respect to the CO (coordinate origin of the GNC system). If r_bp = [0 0 0]' 
% the CG coincides with the CO and the matrices MRB and CRB are computed 
% in the CO = CG. Else if the distance vector from the CO to the CG is 
% r_bg = [xg, yg, zg]' use r_bp = r_bg.
%
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,[0, 0, 0]') computes MRB and CRB in CG
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bg) computes MRB and CRB in CO
%
% In general nu = [u, v, w, p, q, r]' but for a rigid body there exists a
% CRB matrix which only depends on the angular velocity nu2 = [p, q, r,]'. 
% In other words, CRB is independent of linear velocity nu1 = [u, v, w]'. 
% This is particular useful for the relative equations of motion where 
% nu_r = nu - nu_c and nu_c = [u_c, v_c, w_c, 0, 0, 0]' is the irrotational
% current velocity. This implies that the marine craft equations of motion
% expressed in the CO satisfies (Fossen 2021, Chapter 3.3.1):
% 
%  MRB * nudot + CRB(nu) * nu = MRB * nudot_r + CRB(nu_r) * nu_r = tau 
%
% Outputs:
%  MRB                 - Rigid-body mass matrix
%  CRB = CRB(nu2)      - Coriolis-centripetal matrix, independet of nu1
%
% Inputs:
%  nu2 = [p, q, r]'    - Angular velocity vector 
%  R44, R55, R66       - Radii of gyration with respect to CG
%  r_bp                - Vector from CO to an arbitrarily point P, e.g.
%                        r_bp = [0, 0, 0]'    CO = CG
%                        r_bp = r_bg = [xg, yg, zg]' vector from CO to CG
% 
% Author:    Thor I. Fossen
% Date:      2019-04-19
% Revisions: 
%  2021-04-21 - Improved the documentation

I3 = eye(3);
O3 = zeros(3,3);

% MRB and CRB expressed in the CG
I_g = m * diag([R44^2, R55^2, R66^2]);
MRB_CG = [ m * I3  O3
          O3      I_g ];
CRB_CG = [ m * Smtrx(nu2)    O3
           O3               -Smtrx(I_g*nu2)  ];
      
% Transform the matrices MRB and CRB from the CG to the CO           
H = Hmtrx(r_bp);      
MRB = H' * MRB_CG * H;
CRB = H' * CRB_CG * H;

