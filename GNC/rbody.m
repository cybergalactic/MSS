function [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_p)
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_p) computes the 6x6 rigid-body mass
% and Coriolis-centripetal matrices for an arbitrariliy point P. 
%
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,[0, 0, 0]') computes MRB and CRB in CG.
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_g) computes MRB and CRB in CO.
%
%Equations of motion:
%    MRB * nu_dot + CRB(nu2) * nu = tau, where nu = [u,v,w,p,q,r]'
%
% Outputs:
% MRB                 - Rigid-body mass matrix
% CRB(nu2)            - Coriolis-centripetal matrix, independet of nu1
%
% Inputs:
% nu2 = [p, q, r]'    - Angular velocity vector 
% R44, R55, R66       - Radii of gyration with respect to CG
% r_p                 - Vector from CO to an arbitrarily point P, e.g.
%                       r_p = [0, 0, 0]'    CO = CG
%                       r_p = r_g = [xg, yg, zg]' vector from CO to CG
% 
% Author:    Thor I. Fossen
% Date:      19 April 2019
% Revisions: 

I3 = eye(3);
O3 = zeros(3,3);

% MRB and CRB in CG
I_g = m * diag([R44^2, R55^2, R66^2]);
MRB_CG = [ m * I3  O3
          O3      I_g ];
CRB_CG = [ m * Smtrx(nu2)    O3
           O3               -Smtrx(I_g*nu2)  ];
      
% Transform MRB and CRB from CG to CO           
H = Hmtrx(r_p);      
MRB = H' * MRB_CG * H;
CRB = H' * CRB_CG * H;

