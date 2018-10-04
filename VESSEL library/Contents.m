% MSS Vessel Models
%
% Standard ship, semi-submersible and underwater vehicle models.
%
%    Ships:
%    mariner     - Mariner class vessel, L=160 m (nonlinear maneuvering model).
%    tanker      - Esso Osaka tanker, L=304 m (nonlinear course unstable maneuvering model).
%    container   - Container ship, L=175 m (nonlinear maneuvering model including the roll mode)
%    Lcontainer  - Container ship, L=175 m (LINEAR maneuvering model including the roll mode)
%    navalvessel - Multipurpose naval vessel, L = 51.5 m (nonlinear manneuvering model)
%    supply      - Supply vessel, L = 76.2 m (linear DP model)
%
%    Underwater vehicles:
%    DSRV        - Deep submergence rescue vehicle (DSRV), L = 5.0 m
%    npsauv      - Naval Postgraduate School autnomous underwater vehicle (AUV), L = 5.3 m
%
%    Semi-submersibles
%    rig         - Semi-sub MDG (mass-damper-spring) model, L = 84.6 m
%
%    Time-series simulations:
%    SIMmariner.m   - Simulate mariner.m under PD control
%    SIMcontainer.m - SImulate container.m and Lcontainer.m under PD control
%
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org
