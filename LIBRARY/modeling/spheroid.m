function [MRB,CRB] = spheroid(a,b,nu2,r_bg)
% [MRB,CRB] = spheroid(a,b,nu2,r_bg) computes the 6x6 rigid-body mass
% and Coriolis-centripetal matrices of a prolate spheroid of length
% L = 2 * a and diameter D = 2 * b. The spheroid can be used to approximate 
% a cylinder-shaped autonomous underwater vehicle (AUV). In general 
% nu = [u,v,w,p,q,r]', while linear and angular velocities are denoted by
% nu1 = [u, v, w]' and nu2 = [p, q, r]'. The CRB matrix is computed using 
% the linear velocity-independent representation (Fossen 2021, Chapter 3.3.1)
% according to:
%
%  [MRB,CRB] = spheriod(a, b,[p, q, r]',[xg, yg, zg]') 
%
% This is particular useful for the relative equations of motion where 
% nu_r = nu - nu_c and nu_c = [u_c, v_c, w_c, 0, 0, 0]' is the irrotational
% current velocity. This implies that the AUV equations of motion
% expressed in the CO satisfies:
% 
%  MRB * nudot + CRB(nu) * nu = MRB * nudot_r + CRB(nu_r) * nu_r = tau
%
% Outputs:
%  MRB:                 Rigid-body mass matrix
%  CRB = CRB(nu2):      Coriolis-centripetal matrix, independent of nu1
%
% Inputs:
%  a, b:                Semiaxes a > b
%  nu2 = [p, q, r]':    Angular velocity vector 
%  r_bg:                r_bg = [xg, yg, zg]' vector from the CO to the CG
% 
% Author:    Thor I. Fossen
% Date:      24 April 2021 

O3 = zeros(3,3);

% Mass of spheriod 
rho = 1025;
m = 4/3 * pi * rho * a * b^2;   

% Moment of inertia
Ix = (2/5) * m * b^2;
Iy = (1/5) * m * (a^2 + b^2);
Iz = Iy;
Ig = diag([Ix Iy Iz]);

% Rigid-body matrices expressed in the CG
MRB_CG = diag([ m m m Ix Iy Iz ]);
CRB_CG = [ m * Smtrx(nu2)    O3
           O3               -Smtrx(Ig*nu2) ];
       
% Transform MRB and CRB from the CG to the CO        
H = Hmtrx(r_bg);
MRB = H' * MRB_CG * H;    
CRB = H' * CRB_CG * H;


       



      

