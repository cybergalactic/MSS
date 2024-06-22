% exShipHydrostatics is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script calculates various ship parameters including dimensions, 
% physical properties, center of gravity (CG), center of buoyancy (CB), 
% moments of inertia, metacentric heights, and the G matrix based on 
% Fossen (2021 Chapter 4.2).
%
% References: 
%      T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%           Motion Control. 2nd. Edition, Wiley. URL: www.fossen.biz/wiley 
%
% Author:    Thor I. Fossen
% Date:      2024-07-22
% Revisions: 

% Ship main characteristics
L = 92;                     % Length (m)
B = 21;                     % Beam (m)
T = 6;                      % Draft (m)

% Physical constants
rho = 1025;                 % Density of water (kg/m^3)

% Mass properties and block coefficient: Cb = nabla / (L * B * T)
Cb = 0.75;                  % Block coefficient, dimensionless
nabla = Cb * L * B * T;     % Volume displacement (m^3)
m = rho * nabla;            % Mass (kg)

% CG location relative to the midships coordinate origin (CO)
r_bg = [-0.5 0 -1]';        % CG position (m)

% Waterplane area coefficient: Cw = Awp / (L * B)
Cw = 0.85;                  % Waterplane area coefficient, dimensionless
Awp = Cw * B * L;           % Waterplane area (m^2)

% Calculating KB (height of the center of buoyancy above keel)
KB = (1/3) * (5*T/2 - nabla/Awp);    % Equation (4.38)

% Munro-Smith coefficient calculation
k_munro_smith =  (6 * Cw^3) / ((1+Cw) * (1+2*Cw));    % Equation (4.37)

% Center of buoyancy (CB) location relative to the coordinate origin (CO)
r_bb = [-0.5 0 T-KB]';      % CB position (m)

% Vertical distance between CG and CB
BG = r_bb(3) - r_bg(3);     % Vertical distance (m)

% Moments of inertia
I_T = k_munro_smith * (B^3 * L) / 12;   % Transverse moment of inertia (m^4)
I_L = 0.7 * (L^3 * B) / 12;             % Longitudinal moment of inertia (m^4)

% Metacentric heights
BM_T = I_T / nabla;         % Transverse metacentric radius (m)
BM_L = I_L / nabla;         % Longitudinal metacentric radius (m)

% Metacentric heights (should be between 0.5 and 1.5 m for transverse)
GM_T = BM_T - BG;           % Transverse metacentric height (m)
GM_L = BM_L - BG;           % Longitudinal metacentric height (m)

% G matrix calculation
LCF = -0.5;                 % x-distance from the CO to the center of Awp
r_bp = [0 0 0]';            % Zero vector for G matrix computation

% Compute the G matrix in the CO
G = Gmtrx(nabla, Awp, GM_T, GM_L, LCF, r_bp);

%% PRINT SHIP DATA
fprintf('%s\n','-------------------------------------------------------------------------------------');
fprintf('%s\n','SHIP HYDROSTATICS');
fprintf('%s\n','-------------------------------------------------------------------------------------');
fprintf('%-40s %8.2f m \n', 'Length (L):', L);
fprintf('%-40s %8.2f m \n', 'Beam (B):', B);
fprintf('%-40s %8.2f m \n', 'Draft (T):', T);
fprintf('%-40s %8.2f kg \n', 'Mass (m):', m);
fprintf('%-40s %8.2f kg/m^3 \n', 'Density of water (rho):', rho);
fprintf('%-40s %8.2f m^3 \n', 'Volume displacement (nabla):', nabla);
fprintf('%-40s %8.2f \n', 'Block coefficient (C_b):', Cb);
fprintf('%-40s %8.2f \n', 'Waterplane area coefficient (C_w):', Cw);
fprintf('%-40s [%2.1f %2.1f %2.1f] m \n', 'Center of gravity (r_bg):',...
    r_bg(1), r_bg(2), r_bg(3));
fprintf('%-40s %8.2f m \n', 'Transverse metacentric height (GM_T):', GM_T);
fprintf('%-40s %8.2f m \n', 'Longitudinal metacentric height (GM_L):', GM_L);

% Print the G matrix
fprintf('%s\n','-------------------------------------------------------------------------------------');
fprintf('%s\n','G MATRIX');
fprintf('%s\n','-------------------------------------------------------------------------------------');
fprintf('%s\n', 'G = ');
disp(G);
