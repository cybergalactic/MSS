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
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

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

