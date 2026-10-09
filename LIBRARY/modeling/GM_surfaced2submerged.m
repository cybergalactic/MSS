function [GM, BM, zB] = GM_surfaced2submerged( ...
    I_waterplane, nabla, zn, T, zB_surface, zB_submerged, zG)
% Computes the transverse or longitudinal metacentric height (GM), metacentric
% radius (BM), and the center of buoyancy zB for an underwater vehicle
% based on the submersion level, zn (Fossen 2027, Chapter 4).
%
%   BG = -(zG - zB)
%   GM = BM - BG = BM + zG - zB
%   GM = BM + KB - KG                    - Alternative formula
%
%   alpha = exp(-5 * (zn / T)^2)         - Smooth transition parameter
%   BM = alpha * (I_waterplane / nabla)  - alpha between 0 and 1
%
% Surfaced vehicle (alpha = 1):
%   GM = (I_waterplane / nabla) + zG - zB
%
% Submerged vehicle (alpha approaches 0):
%   GM = zG - zB                          - Since BM approaches 0
%      = -BG
%
% The exponential law is a smooth transition model. Exact hydrostatic
% calculations require the displacement volume, waterplane moment of inertia,
% and center of buoyancy to be computed from the immersed hull geometry.
%
% INPUTS:
%   I_waterplane  - Moment of inertia of the waterplane area. This can be either 
%                   the transverse moment of inertia (I_T) or longitudinal moment 
%                   of inertia (I_L), depending on the stability axis.
%   nabla         - Displacement volume of the vehicle.
%   zn            - Current depth of the vehicle below the waterline.
%   T             - Draft of the vehicle at the surface.
%   zB_surface    - Downwards position of the CB (zB) relative to the CO when
%                   the vehicle is at the surface.
%   zB_submerged  - Downwards position of the CB (zB) relative to the CO when
%                   the vehicle is fully submerged.
%   zG            - Downwards position of the CG (zG) relative to the CO.
%
% OUTPUTS:
%   GM            - The metacentric height for the specified stability axis. 
%                   Returns GM_T for transverse stability when I_waterplane = I_T, 
%                   and GM_L for longitudinal  stability when I_waterplane = I_L.
%   BM            - The metacentric radius for the specified stability axis. 
%                   Returns BM_T for transverse stability when I_waterplane = I_T, 
%                   and BM_L for longitudinal tability when I_waterplane = I_L.
%   zB            - The vertical position of the CB (zB) relative to the CO
%                   based on the submersion level zn, interpolating between 
%                   zB_surface and zB_submerged.
% 
% Example calls:
%
%   [GM_T, BM_T, zB] = GM_surfaced2submerged( ...
%      I_T, nabla, zn, T, zB_surface, zB_submerged, zG)
%
%   [GM_L, BM_L, zB] = GM_surfaced2submerged( ...
%      I_L, nabla, zn, T, zB_surface, zB_submerged, zG)
%
%   See 'exPlotGM.m' for a numerical example showing GM_T as a function of depth.
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:    Thor I. Fossen
% Date:      2024-11-08
% Revisions:

rangeCheck(zn, -10, 10000); % The depth should be between -10 and 10000 m
rangeCheck(T, 0, 100); % The draft should be between 0 and 100 m

% Submersion ratio alpha, between 0 and 1
alpha = exp(-5*(zn/T)^2);

% Compute hydrostatic parameter using alpha as transition parameter
BM = alpha * (I_waterplane / nabla);
zB = (1 - alpha) * zB_submerged + alpha * zB_surface;

% Compute GM_T
BG = -(zG - zB);
GM = BM - BG;

end
