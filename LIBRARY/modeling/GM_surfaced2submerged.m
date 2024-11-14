function [GM, BM, z_b] = GM_surfaced2submerged( ...
    I_waterplane, nabla, zn, T, z_b_surface, z_b_submerged, z_g)
% GM_surfaced2submerged is compatible with MATLAB and GNU Octave (www.octave.org).
% Computes the transverse or longitudinal metacentric height (GM), metacentric
% radius (BM), and the center of buoyancy (z_b) for an underwater vehicle 
% based on the submersion level, zn (Fossen 2021, Chapter 4).
%
%   GM = BM + z_g - z_b
%   GM = BM + KB - KG               - Alternative formula
%
%   BM = (1 - alpha) * (I / nabla)  - alpha between 0 and 1
%
% Surfaced vehicle (alpha = 0):
%   GM = I_waterplane / nabla) + z_g - z_b
%
% Submerged vehicle (alpha = 1):
%   GM = z_g - z_b                    - Since I_waterplane = 0
%      = BG_z 
%
% INPUTS:
%   I_waterplane  - Moment of inertia of the waterplane area. This can be either 
%                   the transverse moment of inertia (I_T) or longitudinal moment 
%                   of inertia (I_L), depending on the stability axis.
%   nabla         - Displacement volume of the vehicle.
%   zn            - Current depth of the vehicle below the waterline.
%   T             - Draft of the vehicle at the surface.
%   z_b_surface   - Downwards position of the CG (z_b) relative to the CO when 
%                   the vehicle is at the surface.
%   z_b_submerged - Downwards position of the CB (z_b) relative to the CO when 
%                   the vehicle is fully submerged.
%   z_g           - Downwards position of the CG (z_g) relative to the CO.
%
% OUTPUTS:
%   GM            - The metacentric height for the specified stability axis. 
%                   Returns GM_T for transverse stability when I_waterplane = I_T, 
%                   and GM_L for longitudinal  stability when I_waterplane = I_L.
%   BM            - The metacentric radius for the specified stability axis. 
%                   Returns BM_T for transverse stability when I_waterplane = I_T, 
%                   and BM_L for longitudinal tability when I_waterplane = I_L.
%   z_b           - The vertical position of the CB (z_b) relative to the CO 
%                   based on the submersion level zn, interpolating between 
%                   z_b_surface and z_b_submerged.
% 
% Example calls:
%
%   [GM_T, BM_T, z_b] = GM_surfaced2submerged( ...
%      I_T, nabla, zn, T, z_b_surface, z_b_submerged, z_g)
%
%   [GM_L, BM_L, z_b] = GM_surfaced2submerged( ...
%      I_L, nabla, zn, T, z_b_surface, z_b_submerged, z_g)
%
%   See 'exPlotGM.m' for a numerical example showing GM_T as a function of depth.
%
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-11-08
% Revisions:

rangeCheck(zn, -10, 10000); % The depth should be between -10 and 1000 m
rangeCheck(T, 0, 100); % The draft should be between 0 and 100 m

% Submersion ratio alpha, between 0 and 1
alpha = exp(-5*(zn/T)^2);

% Compute hydrostatic parameter using alpha as transition parameter
BM = alpha * (I_waterplane / nabla);
z_b = (1 - alpha) * z_b_submerged + alpha * z_b_surface;

% Compute GM_T
GM = BM + (z_g - z_b);

end