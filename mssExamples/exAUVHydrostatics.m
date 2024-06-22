% exAUVHydrostatics is compatible with MATLAB and GNU Octave (www.octave.org).
% This script calculates various AUV parameters including dimensions,
% physical properties, center of gravity (CG), center of buoyancy (CB),
% moments of inertia, metacentric heights, and the g vector based on
% Fossen (2021, Chapter 4.1).
%
% References:
%      T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and
%           Motion Control. 2nd Edition, Wiley. URL: www.fossen.biz/wiley
%
% Author:    Thor I. Fossen
% Date:      2024-07-22
% Revisions: 

% AUV main characteristics
L_auv = 1.6;                % Length (m)
D_auv = 0.19;               % Cylinder diameter (m)

% Physical constants
rho = 1025;                 % Density of water (kg/m^3)
g = 9.81;                   % Gravitational acceleration (m/s^2)

% CG location relative to the midships coordinate origin (CO)
r_bg = [0 0 0.02]';         % CG position (m)

% Center of buoyancy (CB) location relative to the coordinate origin (CO)
r_bb = [0 0 0]';            % CB position (m)

% Rigid-body mass, Fossen (2021, Ch. 8.4.2)
a = L_auv / 2;              % Spheroid semi-axes a and b
b = D_auv / 2; 
[MRB, CRB] = spheroid(a, b, [0 0 0]', r_bg);
m = MRB(1,1);
nabla = m / rho;  % Calculate volume displacement
W = m * g;        % Calculate the weight
B = W;            % Buoyancy = Weight

% User input for phi and theta
phi = input('Enter roll angle in degrees: ');
theta = input('Enter pitch angle in degrees: ');
theta = deg2rad(theta);
phi = deg2rad(phi);

% Compute the restoring forces and moments in the CO 
gVect = gvect(W, B, theta, phi, r_bg, r_bb);

%% PRINT SHIP DATA
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%s\n', 'AUV HYDROSTATICS');
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%-40s %8.2f m \n', 'Length (L):', L_auv);
fprintf('%-40s %8.2f m \n', 'Cylinder diameter (d):', D_auv);
fprintf('%-40s %8.2f kg \n', 'Mass (m):', m);
fprintf('%-40s %8.2f kg/m^3 \n', 'Density of water (rho):', rho);
fprintf('%-40s %8.2f m^3 \n', 'Volume displacement (nabla):', nabla);
fprintf('%-40s %8.2f N \n', 'Weight (W):', W);
fprintf('%-40s %8.2f N \n', 'Buoyancy (B):', B);
fprintf('%-40s [%2.1f %2.1f %2.1f] m \n', 'Center of gravity (r_bg):', r_bg(1), r_bg(2), r_bg(3));
fprintf('%-40s [%2.1f %2.1f %2.1f] m \n', 'Center of buoyancy (r_bb):', r_bb(1), r_bb(2), r_bb(3));

% Print the g vector
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%s\n', 'g VECTOR');
fprintf('%s\n', '-------------------------------------------------------------------------------------');
fprintf('%s\n', 'g = ');
disp(gVect);
