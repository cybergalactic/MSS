function C = m2c(M,nu)
% C = m2c(M,nu) computes the Coriolis-centripetal matrix C(nu) from the
% the system inertia matrix M > 0 for varying velocity nu. 
% If M is a 6x6 matrix and nu = [u, v, w, p, q, r]', the output is a 6x6 C matrix
% If M is a 3x3 matrix and nu = [u, v, r]', the output is a 3x3 C matrix

% Examples:   C_RB = m2c(M_RB,nu)
%             C_A  = m2c(M_A, nu)
% 
% See also: [M_RB,C_RB] = rbody(m,R44,R55,R66,nu2,r_p) 
%
% Author:    Thor I. Fossen
% Date:      14 June 2001
% Revisions: 26 June 2002, M21 = M12 is corrected to M12'
%            10 Jan 2004, the computation of C is generalized to a nonsymmetric M > 0
%            22 Oct 2020, generalized to accept 3-DOF surface craft models

M = 0.5 * (M + M');                   % symmetric inertia matrix

if (length(nu) == 6)                  % 6-DOF model
    
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
    
else                                 % 3-DOF model (surge, sway and yaw)
    
    C = [ 0             0            -M(2,2)*nu(2)-M(2,3)*nu(3)
          0             0             M(1,1)*nu(1)
          M(2,2)*nu(2)+M(2,3)*nu(3)  -M(1,1)*nu(1)   0              ];
      
    
end

