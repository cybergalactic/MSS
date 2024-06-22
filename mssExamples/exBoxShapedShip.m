% exBoxShapedShip is compatible with MATLAB and GNU Octave (www.octave.org).
% This script computes the transverse metacentric height GM_T, and the 
% heave and roll periods of a box-shaped ship with the coordinate origin (CO)
% at midships on the centerline using the hydrostatic formulas by Fossen 
% (2021, Chapter 4.2).
%
% References: 
%  T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%      Motion Control. 2nd Edition, Wiley. URL: www.fossen.biz/wiley 
%
% Author:    Thor I. Fossen
% Date:      2024-07-22
% Revisions: 

% Constants
rho = 1025;       % Water density (kg/m^3)
g = 9.82;         % Acceleration due to gravity (m/s^2)

% Box-shaped ship
L = 80;           % Length (m)
B = 20;           % Beam (m)
T = 6;            % Draft (m)
R44 = 0.35 * B;   % Radius of gyration (m)
zb = T / 2;       % Center of buoyancy (CB) w.r.t. the CO (m)
zg = -2;          % Center of gravity (CG) w.r.t. the CO (m)

% Displaced volume and mass
Cb = 1.0;                % Block coefficient, dimensionless
nabla = Cb * L * B * T;  % Displaced volume (m^3)
m = rho * nabla;         % Displacement mass (kg)

% Hydrostatics (Chapter 4.2.3)
I_T = (1/12) * B^3 * L;  % Transverse moment of inertia (m^4)
BM_T = I_T / nabla;      % Transverse metacentric radius (m)
BG = zb - zg;            % Vertical distance between CB and CG (m)
GM_T = BM_T - BG;        % Transverse metacentric height (m)

Ix = m * R44^2;          % Moment of inertia about the CG (kg·m^2)
Ix_CF = Ix + m * zg^2;   % Moment of inertia about the CF (kg·m^2)
A44_CF = 0.2 * Ix_CF;    % Added moment of inertia (kappa4 = 0.2) (kg·m^2)

% Heave and roll periods in seconds w.r.t. the CF (Chapter 4.3.3)
T3 = 2 * pi * sqrt( 2 * T / g );  
T4 = 2 * pi * sqrt( (Ix_CF + A44_CF) / (rho * g * nabla * GM_T) ); 

%% PRINT DATA
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%s\n', 'BOX-SHAPED SHIP');
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%-40s %8.2f m \n', 'Length (L):', L);
fprintf('%-40s %8.2f m \n', 'Beam (B):', B);
fprintf('%-40s %8.2f m \n', 'Draft (T):', T);
fprintf('%-40s %8.2f kg \n', 'Mass (m):', m);
fprintf('%-40s %8.2f kg/m^3 \n', 'Density of water (rho):', rho);
fprintf('%-40s %8.2f m^3 \n', 'Volume displacement (nabla):', nabla);
fprintf('%-40s %8.2f \n', 'Block coefficient (C_b):', Cb);
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%s\n', 'HYDROSTATICS');
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%-40s %8.2f m \n', 'Transverse metacentric height (GM_T):', GM_T);
fprintf('%-40s %8.2f s \n', 'Heave period (T3):', T3);
fprintf('%-40s %8.2f s \n', 'Roll period (T4):', T4);