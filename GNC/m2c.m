function C = m2c(M,nu)
% C = m2c(M,nu) computes the 6x6 Coriolis-centripetal matrix C(nu) from the
% the 6x6 system inertia matrix M > 0 using the velocity vector
% nu = [u, v, w, p, q, r]'.
%
% Examples:   C_RB = m2c(M_RB,nu)
%             C_A  = m2c(M_A, nu)
% 
% See also: [M_RB,C_RB] = rb(m,RX,RY,RZ,nu2,r_p) 
%
% Author:    Thor I. Fossen
% Date:      14 June 2001
% Revisions: 26 June 2002, M21 = M12 is corrected to M12'
%            10 Jan 2004, the computation of C is generalized to a nonsymmetric M > 0

% Kinetic energy: 
%    T = 1/2 (v'Mv) = 1/2 v'(1/2(M + M')+1/2(M-M')) v =  1/2 v' Msym v
%    Msym = 1/2 (M+M')
Msym = 0.5*(M+M');
M = Msym;

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

