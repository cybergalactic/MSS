% exPlotGM
% Example script to compute and plot hydrostatic stability parameters (GM, BM)
% and center of buoyancy (z_b) for an underwater vehicle transitioning from
% surfaced to fully submerged conditions.
%
% [GM_T, BM_T, z_b] = GM_surfaced2submerged( I_T, nabla, zn, T, z_b_surface, ... 
%    z_b_submerged, z_g) computes the transverse stability parameters.
%
% Author:    Thor I. Fossen
% Date:      2024-11-08
% Revisions:
close all

% Initialize vessel parameters and properties
targetDepth = 10;               % Maximum plotting depth
L = 4;                          % Length of the vessel
B = 2;                          % Beam of the vessel
T = 2;                          % Draft of the vessel
cB = 0.6;                       % Block coefficient
I_T = 1/12 * B^3 * L;           % Transverse moment of inertia of the waterplane area
nabla = cB * L * B * T;         % Displacement volume
z_b_surface = (1/3) * T;        % Center of buoyancy at surface
z_b_submerged = T / 2;          % Center of buoyancy when fully submerged
z_g = T / 2 + 0.3;              % Center of gravity

% Initialize arrays to store values for plotting
zn_values = 0:0.2:targetDepth;
BG_z_values = zeros(size(zn_values));
GM_T_values = zeros(size(zn_values));
BM_T_values = zeros(size(zn_values));
z_g_values = z_g * ones(size(zn_values)); % z_g is constant
z_b_values = zeros(size(zn_values));

% Populate arrays with computed values
N = length(zn_values);
for i = 1:N
    zn = zn_values(i);
    [GM_T, BM_T, z_b] = GM_surfaced2submerged( ...
        I_T, nabla, zn, T, z_b_surface, z_b_submerged, z_g);
    BG_z_values(i) = z_g - z_b;
    GM_T_values(i) = GM_T;
    BM_T_values(i) = BM_T;
    z_b_values(i) = z_b;
end

% Plot solid lines for each parameter
plot(zn_values, BG_z_values, 'g-', ...
    zn_values, GM_T_values, 'r-', ...
    zn_values, BM_T_values, 'b-')
hold on
plot(zn_values, z_g_values, 'co-', 'MarkerIndices', 1:5:N)
plot(zn_values, z_b_values, 'ko-', 'MarkerIndices', 1:5:N)
hold off
grid on

legend('BG_z = z_g - z_b', 'GM_T = BM_T + BG_z', ...
    'BM_T = exp(-(zn/T)^2) * (I_T / nabla)', ...
    'z_g (Positive downwards)', 'z_b (Positive downwards)', ...
    'location','E')
title('Hydrostatic Stability Parameters as a Function of Submersion Depth')
xlabel('Depth Transition: Surfaced to Submerged')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',12)