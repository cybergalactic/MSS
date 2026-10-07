function C = m2c(M,nu)
% C = m2c(M,nu) computes the Coriolis-centripetal matrix C(nu) from the
% the system inertia matrix M > 0 for varying velocity nu. 
% If M is a 6x6 matrix and nu = [u, v, w, p, q, r]', the output is a 6x6 C matrix
% If M is a 3x3 matrix and nu = [u, v, r]', the output is a 3x3 C matrix.
%
% Examples: CRB = m2c(MRB, nu)     
%           CA  = m2c(MA, nu)
% Output:
%  C:  Coriolis-centripetal matrix C = C(nu) 
%
% Inputs:
%  M:  6x6 or 3x3 rigid-body MRB or added mass MA system matrix 
%  nu: nu = [u, v, w, p, q, r]' or nu = [u, v, r]'
%
% The Coriolis and centripetal matrix depends on nu1 = [u,v,w]' and nu2 =
% [p,q,r]' as shown in Fossen (2027, Theorem 3.2). Alternatively, the matrix
% CRB = CRB(nu2) can be computed using 
% 
% [MRB, CRB] = rbody(m, R44, R55, R66, nu2, r_bP) 
%
% which only requires the angular velocity vector nu2 = [p, q, r]'. This is 
% known as the linear velocity-independent representation.
%
% Author:    Thor I. Fossen
% Date:      2001-06-14
% Revisions: 2002-06-26  M21 = M12 is corrected to M12'.
%            2004-01-10  The computation of C(nu) is generalized to a 
%                        nonsymmetric M > 0 (experimental data).
%            2020-10-22  Generalized to accept 3-DOF horizontal-plane models.
%            2021-04-24  Updated the documentation.
%            2026-10-07  Corrected the 3-DOF formulation to use the full
%                        momenta (E. Krizamn)

M = 0.5 * (M + M');      % Symmetrization of the inertia matrix

if (length(nu) == 6)     % 6-DOF model
     
    M11 = M(1:3,1:3);
    M12 = M(1:3,4:6);
    M21 = M12';
    M22 = M(4:6,4:6);

    nu1 = nu(1:3);
    nu2 = nu(4:6);
    nu1_dot = M11 * nu1 + M12 * nu2;
    nu2_dot = M21 * nu1 + M22 * nu2;

    C = [  zeros(3,3)      -Smtrx(nu1_dot)
          -Smtrx(nu1_dot)  -Smtrx(nu2_dot) ];
    
else   % 3-DOF model (surge, sway and yaw)
    
    p = M * nu;
    C = [ 0     0    -p(2)
          0     0     p(1)
          p(2) -p(1)  0    ];
    
end

