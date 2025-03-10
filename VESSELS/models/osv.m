function [xdot,U,M] = osv(x,ui,Vc,betaVc)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org)
% [xdot,U,M] = osv(x,ui,Vc,betaVc) returns the speed U in m/s (optionally) 
% and the time derivative xdot of the state vector for an Offshore Supply 
% vessel (OSV). The 6x6 mass matrix M is an optionally output, which can be 
% used for control design.The 6-DOF equations of motion arebased on the 
% nonlinear model of Fossen (2021, Eqs. 6.111-6.116) given by
%   
%   eta_dot = J(eta) * nu
%   nu_dot = nu_c_dot + Minv * ( tau_thr +  tau_drag + tau_crossflow...
%          - (CRB + CA + D) * nu_r - G * eta )
%
% where Minv = inv(MRB + MA) and nu_r = nu - nu_c is the relative 
% velocity. The OSV equipped by two bow tunnel thrusters and two stern 
% azimuth thrusters. 
%
% Inputs: 
%   x: state vector x = [ u v w p q r x y z phi theta psi ]'
%   ui: control inputs ui = [n1, n2, n3, n4, alpha1, alpha2] where n1, n2, 
%       n3 and n4 are the propeller speeds (rps). The last two arguments 
%       alpha1 and alpha2 are the azimuth angles (rad)
%   Vc: OPTIONAL current speed
%   betaVc: OPTIONAL current direction (rad)
%
%   The arguments Vc (m/s) and betaVc (rad) are optional arguments for 
%   ocean currents given by
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), 0 ] 
% 
% The generalized thrust vector satisfy (Fossen 2021, Section 11.2.1)
%
%   tau_thr = T_thr(alpha) * K_thr * u_thr
%
% where
%
%   u_thr = [bowThruster 1, bowThruster 2, sternAzimuth 1, sternAzimuth 2]
%   alpha = [alpha1, alpha2]'    
%
% and u_thr(i) = abs(ui(i)) * ui(i) for i = 1...4 is the squared propeller 
% speed. The azimuth angles are defined by alpha = [ui(5), ui(6)]'. The OSV 
% main characteristics are displayed by calling the function OSV without 
% input arguments. The function calls are:
%
%   osv();                                 : Display the OSV main data
%   [~,~,M] = osv()                        : Return the 6x6 mass matrix M
%   [xdot,U] = osv(x,ui,Vc,betaVc,alphaVc) : Return xdot and U, 2-D currents
%   [xdot,U] = osv(x,ui)                   : Return xdot and U, no currents 
% 
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:
%   2024-04-22: Enhanced compatibility with GNU Octave
%   2024-06-07: Added M as an optional output argument

persistent vessel;

% Initialization of vessel parameters. Compute the vessel parameters only once
% to avoid that RK4 repeats the computation in the loop at each time step
if isempty(vessel)

    %% Ship model parameters
    vessel.L = 83;                 % Length (m)
    vessel.B = 18;                 % Beam (m)
    vessel.T = 5;                  % Draft (m)
    vessel.rho = 1025;             % Density of water (kg/m3)
    vessel.Cb = 0.65;              % Block coefficient: Cb = nabla / (L * B * T)
    vessel.S = vessel.L * vessel.B + 2 * vessel.T * vessel.B; % Wetted surface, box approximation

    % Thrust: T = K_max * abs(n/n_max) * (n/n_max) = K_thr * abs(n) * n
    vessel.K_max = [300e3 300e3 420e3 655e3]';% Max propeller thrust (N)
    vessel.n_max = [140 140 150 200]'; % Max propeller speed (rpm)
    vessel.K_thr = diag(vessel.K_max./vessel.n_max.^2);% Thruster coefficient matrix
    vessel.l_x = [37, 35, -vessel.L/2, -vessel.L/2]; % Thruster x-coordinates
    vessel.l_y = [0, 0, 7, -7]; % Thruster y-coordinates

    vessel.thrust_max = vessel.K_max(3)+vessel.K_max(4); % max thrust in the surge direction (N)
    vessel.U_max = 7.7; % Max cruise speed (m/s) corresponding to max thrust

    vessel.nabla = vessel.Cb * vessel.L * vessel.B * vessel.T; % Volume displacement(m3)
    vessel.m = vessel.rho * vessel.nabla; % Mass (kg)
    vessel.r_bg = [-4.5 0 -1.2]'; % Location of the CG with respect to the CO

    vessel.Cw = 0.8; % Waterplane area coefficient: Cw = Awp/(L * B)
    vessel.Awp = vessel.Cw * vessel.B * vessel.L; % Waterplane area
    vessel.KB = (1/3) * (5*vessel.T/2 - vessel.nabla/vessel.Awp); % Eq. (4.38)
    vessel.k_munro_smith = (6 * vessel.Cw^3) / ...
        ( (1+vessel.Cw) * (1+2*vessel.Cw)); % Eq. (4.37)
    vessel.r_bb = [-4.5 0 vessel.T-vessel.KB]'; % Location of the CB with respect to the CO
    vessel.BG = vessel.r_bb(3) - vessel.r_bg(3); % Vertical distance between CG and CB

    vessel.I_T = vessel.k_munro_smith * ...
        (vessel.B^3 * vessel.L) / 12; % Transverse moment of inertia
    vessel.I_L = 0.7 * (vessel.L^3 * vessel.B) / 12;% Longitudinal moment of inertia
    vessel.BM_T = vessel.I_T / vessel.nabla;
    vessel.BM_L = vessel.I_L / vessel.nabla;
    vessel.GM_T = vessel.BM_T - vessel.BG; % Should be larger than 0.5 m
    vessel.GM_L = vessel.BM_L - vessel.BG;

    % G matrix
    vessel.LCF = -0.5;% x-distance from the CO to the center of Awp
    vessel.r_bp = [0 0 0]'; % Compute G in the CO
    vessel.G = Gmtrx(vessel.nabla,vessel.Awp,vessel.GM_T,vessel.GM_L,vessel.LCF,vessel.r_bp);

    % Rigid-body mass matrix MRB
    vessel.R44 = 0.35 * vessel.B; % Radius of gyration in roll, see Eq.(4.77)-(4.78)
    vessel.R55 = 0.25 * vessel.L; % Radius of gyration in pitch
    vessel.R66 = 0.25 * vessel.L; % Radius of gyration in yaw
    vessel.MRB = rbody(vessel.m,vessel.R44,vessel.R55,vessel.R66,[0,0,0]',vessel.r_bg); 

    % The added mass matrix MA is derived from the frequency-dependent potential
    % coefficients using a look-alike supply vessel in the MSS toolbox. The data
    % is stored in the structure <vessel>.
    %   load supply             % Check data by typing vessel
    %   disp(vessel.main)
    %   vessel = computeManeuveringModel(vessel, 1, 7, [3, 1.2, 3.3], 1);
    vessel.MA = 1e9 * [0.0006   0    0    0    0    0
         0    0.0020         0    0.0031         0   -0.0091
         0         0    0.0083         0    0.0907         0
         0    0.0031         0    0.0748         0   -0.1127
         0         0    0.0907         0    3.9875         0
         0   -0.0091         0   -0.1127         0    1.2416];

    % Mass matrix including hydrodynamic added mass
    vessel.M = vessel.MRB + vessel.MA;
    vessel.Minv = invQR(vessel.M);

    % D matrix
    vessel.T1 = 100;               % Time constants for linear damping (s)
    vessel.T2 = 100;
    vessel.T6 = 1;
    vessel.zeta4 = 0.15;           % Relative damping ratio in roll
    vessel.zeta5 = 0.3;            % Relative damping ratio in pitch
    vessel.D = Dmtrx([vessel.T1, vessel.T2, vessel.T6], ...
        [vessel.zeta4,vessel.zeta5],vessel.MRB,vessel.MA,vessel.G);

end

% Flag for plotting of the surge resitance, linear + quadratic damping
if nargin > 0
    flag = 0;                              
else
    flag = 1;
end

% Initialize state variables
nu = zeros(6,1);
eta = zeros(6,1);

if nargin == 2 || nargin == 0
    Vc = 0;
    betaVc = 0;
end

if nargin == 0                          % Display main ship characteristics
    ui = zeros(6,1);
else
    nu = x(1:6);                        % Generalized velocity vector
    eta = x(7:12);                      % Generalized position vector
end

% Output arguments
M = vessel.M; % Mass matrix
U = sqrt( nu(1)^2 + nu(2)^2 ); % Speed

% Ocean current
v_c = [ Vc * cos(betaVc - eta(6))
        Vc * sin(betaVc - eta(6))
        0                          ];
nu_c = [v_c' zeros(1,3) ]';
nu_c_dot = [-Smtrx(nu(4:6)) * v_c
             zeros(3,1)           ];
nu_r = nu - nu_c;

% Coriolis matrices
[~,CRB] = rbody(vessel.m,vessel.R44,vessel.R55,vessel.R66,nu(4:6),vessel.r_bg'); 
CA  = m2c(vessel.MA, nu);   

% Add linear and quadratic drag in surge using the blending function sigma
[X,Xuu,Xu] = forceSurgeDamping(flag,nu_r(1),vessel.m,vessel.S,vessel.L, ...
    vessel.T1,vessel.rho,vessel.U_max,vessel.thrust_max);
tau_drag = [ X; zeros(5,1) ];

% Avoid double counting, linear and quadratic damping terms
vessel.D(1,1) = 0; % using: X = sigma * Xu * u_r + (1 - sigma) * Xuu * abs(u_r)*u_r

% Add crossflow drag
tau_crossflow = crossFlowDrag(vessel.L,vessel.B,vessel.T,nu_r);

% Thrust: T = K_max * u, where u = abs(n/n_max) * (n/n_max) is normalized
u_thr = abs(ui(1:4)) .* ui(1:4);               % Quadratic propeller speed
alpha = ui(5:6);                               % Azimuth angles
T_thr = thrConfig( {'T', 'T', alpha(1), alpha(2)}, vessel.l_x, vessel.l_y);
tau_3dof = T_thr * vessel.K_thr * u_thr;
tau_thr = [ tau_3dof(1) tau_3dof(2) 0 0 0 tau_3dof(3) ]';

% Kinematics
J = eulerang(eta(4),eta(5),eta(6));

% Equations of motion (Fossen 2021, Eqs. 6.111-6.116)
eta_dot = J * nu;
nu_dot = nu_c_dot + ...
    vessel.Minv * ( tau_thr +  tau_drag + tau_crossflow...
    - (CRB + CA + vessel.D) * nu_r - vessel.G * eta);

xdot = [nu_dot; eta_dot];    

%% Print vessel data
if nargin == 0 && nargout == 0

    % Natural frequencies
    w3 = sqrt( vessel.G(3,3) / vessel.M(3,3) );
    w4 = sqrt( vessel.G(4,4) / vessel.M(4,4) );
    w5 = sqrt( vessel.G(5,5) / vessel.M(5,5) );
    T3 = 2 * pi / w3;
    T4 = 2 * pi / w4;    
    T5 = 2 * pi / w5;

    fprintf('\n');
    fprintf('%s\n','-------------------------------------------------------------------------------------');
    fprintf('%s\n','OFFSHORE SUPPLY VESSEL MAIN CHARACTERISTICS');
    fprintf('%s\n','-------------------------------------------------------------------------------------');
    fprintf('%-40s %8.2f m \n', 'Length (L):', vessel.L);
    fprintf('%-40s %8.2f m \n', 'Beam (B):', vessel.B);
    fprintf('%-40s %8.2f m \n', 'Draft (T):', vessel.T);
    fprintf('%-40s %8.2f kg \n', 'Mass (m):', vessel.m);
    fprintf('%-40s %8.2f kg/m^3 \n', 'Density of water (rho):', vessel.rho);
    fprintf('%-40s %8.2f m^3 \n', 'Volume displacement (nabla):', vessel.nabla);
    fprintf('%-40s %8.2f \n', 'Block coefficient (C_b):', vessel.Cb);
    fprintf('%-40s %8.2f \n', 'Waterplane area coefficient (C_w):', vessel.Cw);
    fprintf('%-40s [%2.1f %2.1f %2.1f] m \n', 'Center of gravity (r_bg):',...
        vessel.r_bg(1), vessel.r_bg(2), vessel.r_bg(3));
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in roll:', vessel.zeta4);
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in pitch:', vessel.zeta5);
    fprintf('%-40s %8.2f m \n', 'Transverse metacentric height (GM_T):', vessel.GM_T);
    fprintf('%-40s %8.2f m \n', 'Longitudinal metacentric height (GM_L):', vessel.GM_L);
    fprintf('%-40s %8.2f s \n', 'Natural period in heave (T3):', T3);
    fprintf('%-40s %8.2f s \n', 'Natural period in roll (T4):', T4);
    fprintf('%-40s %8.2f s \n', 'Natural period in pitch (T5):', T5);   
    fprintf('%-40s %8.2f \n', 'Linear surge damping coefficient (Xu):', Xu); 
    fprintf('%-40s %8.2f \n', 'Quadratic drag coefficient (X|u|u):', Xuu); 
   
    vessel.D(1,1) = -Xu;
    matrices = {'Mass matrix: M = MRB + MA', vessel.M;...
        'Linear damping matrix: D', vessel.D; 'Restoring matrix: G', vessel.G};
    
    for k = 1:size(matrices, 1)
        fprintf('%s\n','-------------------------------------------------------------------------------------');
        fprintf('%-40s\n', matrices{k, 1});
        for i = 1:size(matrices{k, 2}, 1)
            for j = 1:size(matrices{k, 2}, 2)
                if matrices{k, 2}(i,j) == 0
                    fprintf('         0 ');
                else
                    fprintf('%10.2e ', matrices{k, 2}(i,j));
                end
            end
            fprintf('\n');
        end
    end

end

end