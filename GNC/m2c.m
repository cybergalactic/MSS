function C = m2c(M,nu)
% C = m2c(M,nu) computes the Coriolis-centripetal matrix C(nu) from the
% the system inertia matrix M > 0 for varying velocity nu. 
% If M is a 6x6 matrix and nu = [u, v, w, p, q, r]', the output is a 6x6 C matrix
% If M is a 3x3 matrix and nu = [u, v, r]', the output is a 3x3 C matrix.
%
% Examples: CRB = m2c(MRB,nu)     
%           CA  = m2c(MA, nu)
% Output:
%  C         - Coriolis-centripetal matrix C = C(nu) 
%
% Inputs:
%  M        - 6x6 or 3x3 rigid-body MRB or added mass MA system marix 
%  nu       - nu = [u, v, w, p, q, r]' or nu = [u, v, r]'
%
% The Coriolis and centripetal matrix depends on nu1 = [u,v,w]' and nu2 =
% [p,q,r]' as shown in Fossen (2021, Theorem 3.2). It is possible to
% compute C = C(nu2) where nu2 = [p,q,r]' using the linear velocity-
% independent representation, see
%
% [MRB,CRB] = rbody(m,R44,R55,R66,nu2,r_bp) 
%
% Author:    Thor I. Fossen
% Date:      14 Jun 2001
% Revisions: 26 Jun 2002, M21 = M12 is corrected to M12'
%            10 Jan 2004, the computation of C(nU) is generalized to a 
%                         nonsymmetric M > 0 (experimental data)
%            22 Oct 2020, generalized to accept 3-DOF hirizontal-plane models
%            24 Apr 2021, improved the documentation

M = 0.5 * (M + M');      % symmetrization of the inertia matrix

if (length(nu) == 6)     % 6-DOF model
    
    M11 = M(1:3,1:3);
    M12 = M(1:3,4:6);
    M21 = M12';
    M22 = M(4:6,4:6);
    
    nu1 = nu(1:3);
    nu2 = nu(4:6);
    dt_dnu1 = M11*nu1 + M12*nu2;
    dt_dnu2 = M21*nu1 + M22*nu2;
    
    C = [  zeros(3,3)      -Smtrx(dt_dnu1)
        -Smtrx(dt_dnu1)  -Smtrx(dt_dnu2) ];
    
else   % 3-DOF model (surge, sway and yaw)
    
    C = [ 0             0            -M(2,2)*nu(2)-M(2,3)*nu(3)
        0             0             M(1,1)*nu(1)
        M(2,2)*nu(2)+M(2,3)*nu(3)  -M(1,1)*nu(1)   0              ];
    
end

