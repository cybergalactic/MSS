function [M,N] = clarke83(U,L,B,T,Cb,R66,xg,T_surge)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [M,N] = clarke83(U,L,B,T,Cb,R66,xg,T_surge) computes the system matrices 
% of a linear maneuvering model based on Clarke et al. (1983). The  
% hydrodynamic derivatives are based on multiple  linear regression from two 
% sets of model tests. The first data set (Yv, Yr, Nv, Nr) is obtained from 
% rotating arm model experiments, while the second data set 
% (Yvdot, Yrdot, Nvdot, Nrdot, Yv, Yr, Nv, Nr) was obtained from a PMM model.
% Added mass in surge is approximated by Sodings formula: 
%
%   Xudot = -addedMassSurge(m,L)
%
% The time constant in surge is optionally with default value T_surge = L 
% such that:
%
%   Xu = -(m - Xudot) / T_surge
%
% Outputs: 3x3 model matrices M and N in surge, sway and yaw
%      .
%    M nu + N(U) nu = tau,     where N(U) = C(U) + D
% 
% corresponding to the linear maneuvering model
% 
%  (m - Xudot) udot - Xu u                            = (1-t) T
%  (m - Yvdot) vdot + (m - Yrdot)  rdot - Yv v - Yr r = Yd delta
%  (m - Yvdot) vdot + (Iz - Nrdot) rdot - Nv v - Nr r = Nd delta
%
% Note that the coefficients Yv, Yr, Nv and Nr in the N(U) matrix includes 
% linear damping D and the linearized Coriolis and centripetal matrix C(U).
%
% Inputs:
%
%  U:   speed (m/s)
%  L:   length (m)
%  B:   beam (m)
%  T:   draft (m)
%  Cb:  block coefficient (-), Cb = V / (L*B*T) where V is the displaced volume
%  R66: radius of gyration in yaw (smaller vessels R66 ≈ 0.25L, tankers R66 ≈ 0.27L)
%  xg:  x-coordinate of the CG
%  T_surge: (optionally) time constant in surge (default: T_surge = L)
%
% Reference:
%   D. Clarke, P. Gedling, and G. Hine (1983). The application of manoeuvring
%   criteria in hull design using linear theory. Transactions of the Royal
%   Institution of Naval Architects, Vol. 125, pp. 45-68.
% 
% Author:    Thor I. Fossen
% Date:      22 Oct 2020
% Revisions: 
%   2021-06-14 : Removed the C matrix and introduced N(U).
%   2021-12-17 : Xudot is computed by Xudot = -addedMassSurge(m,L).
%   2024-04-19 : Added compability to GNU Octave.

% Rigid-body parameters
rho = 1025;                     % density of water
V = Cb * L * B * T;             % volume displacment
m = rho * V;                    % mass
Iz = m * R66^2 + m * xg^2;      % moment of inerta about the CO

MRB = [ m   0       0           % rigid-body inertia matrix
        0   m       m*xg
        0   m*xg    Iz      ];

% Nondimenisonal hydrodynamic derivatives in surge
Xudot = -addedMassSurge(m,L);
if (nargin == 7)
    T_surge = L; 
end
U = U + 0.001; % avoid singularity for U = 0;
Xu = -((m-Xudot)/T_surge) / (0.5 * rho * L^2 * U);  
Xudot = Xudot / (0.5 * rho * L^3);

% Nondimenisonal hydrodynamic derivatives in sway and yaw
% from Clarke et al. (1983)
S = pi * (T/L)^2;                 % scale factor

Yvdot = -S * ( 1 + 0.16 * Cb * B/T - 5.1 * (B/L)^2 );
Yrdot = -S * ( 0.67 * B/L - 0.0033 * (B/T)^2 );
Nvdot = -S * ( 1.1 * B/L - 0.041 * (B/T) );
Nrdot = -S * ( 1/12 + 0.017 * Cb * (B/T) - 0.33 * (B/L) );
Yv = -S * ( 1 + 0.4 * Cb * (B/T) );
Yr = -S * ( -1/2 + 2.2 * (B/L) - 0.08 * (B/T) );
Nv = -S * ( 1/2 + 2.4 * (T/L) );
Nr = -S * ( 1/4 + 0.039 * (B/T) - 0.56 * (B/L) );

% Nondimenisonal hydrodynamic matrices 
MA_prime = [ -Xudot   0        0
              0      -Yvdot   -Yrdot
              0      -Nvdot   -Nrdot ];

N_prime = [ -Xu   0  0
             0  -Yv -Yr
             0  -Nv -Nr ];
 
% Dimensional model (Fossen 2021, Appendix D)   
T    = diag([1 1 1/L]);
Tinv = diag([1 1 L]);

MA = (0.5 * rho * L^3) * Tinv^2 * (T * MA_prime * Tinv);
N =  (0.5 * rho * L^2 * U) * Tinv^2 * (T * N_prime * Tinv);
 
M = MRB + MA;       % system inertia matrix
