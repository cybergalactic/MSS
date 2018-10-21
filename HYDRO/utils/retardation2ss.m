function  [A,B,C,D] = retardation2ss(K,dT,order)
% RETARDATION2SS (MSS Hydro)
%
% Sys = retardation2ss(K,Ts,order,PlotFlag) estimates a state-space realization 
% from a single DOF retardation function via realization theory.
% The functions uses the Kung's SVD algorithm implementesed in
% imp2ss.m of the robust control toolbox followed by order reduction
% via Square-root balanced truncation balmr.m
%
% Inputs:
%  K     - Samples of the retardation function
%  dT    - Samping period [s]
%  order - Order of the LTI system approximation.
%
% Output:
%  [A,B,C,D] - State-space model
%
%  Author: Tristan Perez
%  Date:      2007-08-13
%  Revisions: 2007-08-24 modified to match Hydro syntax
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

% Identification
scale = max(K);
K = K/scale;
[Ah,Bh,Ch,Dh,TOTBND,SVH] = imp2ss(K,dT,1,1);

% Order reduction
[Ahr,Bhr,Chr,Dhr,totbnd,svh] = balmr(Ah,Bh,Ch,Dh,1,order);

% D2C  
A = Ahr;
B = Bhr;
C = dT*(Chr*scale);
D = dT*(Dhr*scale);

