% exPlotGM
% Example script to compute and plot hydrostatic stability parameters (GM, BM)
% and center of buoyancy zB for an underwater vehicle transitioning from
% surfaced to fully submerged conditions.
%
% [GM_T, BM_T, zB] = GM_surfaced2submerged(I_T, nabla, zn, T, zB_surface, ...
%    zB_submerged, zG) computes the transverse stability parameters.
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
zB_surface = (1/3) * T;         % Center of buoyancy at surface
zB_submerged = T / 2;           % Center of buoyancy when fully submerged
zG = T / 2 + 0.3;               % Center of gravity

% Initialize arrays to store values for plotting
zn_values = 0:0.2:targetDepth;
BG_values = zeros(size(zn_values));
GM_T_values = zeros(size(zn_values));
BM_T_values = zeros(size(zn_values));
zG_values = zG * ones(size(zn_values));  % zG is constant
zB_values = zeros(size(zn_values));

% Populate arrays with computed values
N = length(zn_values);
for i = 1:N
    zn = zn_values(i);
    [GM_T, BM_T, zB] = GM_surfaced2submerged( ...
        I_T, nabla, zn, T, zB_surface, zB_submerged, zG);
    BG_values(i) = -(zG_values(i) - zB);
    GM_T_values(i) = GM_T;
    BM_T_values(i) = BM_T;
    zB_values(i) = zB;
end

% Plot colored lines with distinct markers for grayscale reproduction
markerIndices = 1:5:N;
plot(zn_values, BG_values, 'g-s', ...
    'MarkerIndices', markerIndices, 'MarkerSize', 6)
hold on
plot(zn_values, GM_T_values, 'r-o', ...
    'MarkerIndices', markerIndices, 'MarkerSize', 6)
plot(zn_values, BM_T_values, 'b-^', ...
    'MarkerIndices', markerIndices, 'MarkerSize', 6)
plot(zn_values, zG_values, 'c-d', ...
    'MarkerIndices', markerIndices, 'MarkerSize', 6)
plot(zn_values, zB_values, 'k-x', ...
    'MarkerIndices', markerIndices, 'MarkerSize', 6)
hold off
grid on

legend('BG = -(z_G - z_B)', 'GM_T = BM_T - BG', ...
    'BM_T = exp(-5*(z^n/T)^2) * (I_T / \nabla)', ...
    'z_G (Positive downwards)', 'z_B (Positive downwards)', ...
    'Location', 'best')
title('Hydrostatic Stability Parameters as a Function of Submersion Depth')
xlabel('NED depth z^n (m), positive downwards')

% Fully submerged limiting condition
text(0.65 * targetDepth, ...
    0.5 * (BM_T_values(end) + GM_T_values(end)), ...
    {'Fully submerged:', 'BM_T \approx 0,  GM_T \approx -BG'}, ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'center', ...
    'FontSize', 11)

% Sign conventions: the NED depth coordinate z^n is positive down, while the
% classical hydrostatic distances BG, BM_T, and GM_T are positive up.
ax = gca;
axPosition = get(ax, 'Position');
arrowYLow = axPosition(2) + 0.08;
arrowYHigh = axPosition(2) + 0.24;
arrowXUp = axPosition(1) - 0.075;
arrowXDown = axPosition(1) - 0.04;

% Classical hydrostatic distances: positive upwards
annotation(gcf, 'arrow', [arrowXUp arrowXUp], ...
    [arrowYLow arrowYHigh], ...
    'LineWidth', 1.5, 'HeadLength', 8, 'HeadWidth', 8)
annotation(gcf, 'textbox', ...
    [arrowXUp - 0.055, arrowYHigh + 0.005, 0.11, 0.12], ...
    'String', '+BG  +BM_T +GM_T', 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 12)

% NED depth coordinate: positive downwards
annotation(gcf, 'arrow', [arrowXDown arrowXDown], ...
    [arrowYHigh arrowYLow], ...
    'LineWidth', 1.5, 'HeadLength', 8, 'HeadWidth', 8)
annotation(gcf, 'textbox', ...
    [arrowXDown - 0.05, arrowYLow - 0.05, 0.11, 0.05], ...
    'String', '+z^n', 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 12)

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',11)
