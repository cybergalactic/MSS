function CY_2D = Hoerner(B,T)
% Hoerner computes the 2D Hoerner cross-flow form coefficient as a
% function of B and T. The data is digizied and interpolation is used to
% compute other data points than those in the table
%
%  >> CY_2D = Hoerner(B,T)
%
% Outputs: 
%    CY_2D:    2D Hoerner cross-flow form coefficient
%
% Inputs:
%    T:      draft (m)
%    B:      beam (m)
%
% Author: Thor I. Fossen
% Date:   2007-12-01
%
% Reference:
% A. J. P. Leite, J. A. P. Aranha, C. Umeda and M. B. conti (1998). 
% Current force in tankers and bifurcation of equilibrium of turret
% systems: hydrodynamic model and experiments. 
% Applied Ocean Research 20, pp. 145-256.
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

% CD_DATA = [B/2T  C_D]
CD_DATA = [...    
0.0108623 1.96608
0.176606 1.96573
0.353025 1.89756
0.451863 1.78718
0.472838 1.58374
0.492877 1.27862
0.493252 1.21082
0.558473 1.08356
0.646401 0.998631
0.833589 0.87959
0.988002 0.828415
1.30807 0.759941
1.63918 0.691442
1.85998 0.657076
2.31288 0.630693
2.59998 0.596186
3.00877 0.586846
3.45075 0.585909
3.7379 0.559877
4.00309 0.559315];

CY_2D = interp1(CD_DATA(:,1),CD_DATA(:,2),B/(2*T));
