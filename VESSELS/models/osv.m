function [xdot,U,M] = osv(x,ui,Vc,betaVc)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org)
% [xdot,U,M] = osv(x,ui,Vc,betaVc) returns the speed U in m/s (optionally) 
% and the time derivative xdot of the state vector: 
%
%    x = [ u v w p q r x y z phi theta psi ]' 
%
% for an Offshore Supply vessel (OSV). The 6x6 mass matrix M is an 
% optionally output, which can be used for control design.The 6-DOF 
% equations of motion arebased on the nonlinear model of Fossen (2021, 
% Eqs. 6.111-6.116) given by
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
%   2024-04-22: Enhanced compatibility with GNU Octave.
%   2024-06-07: Added M as an optional output argument

% Flag for plotting of the surge resitance, linear + quadratic damping
if nargin > 0
    flag = 0;                              
else
    flag = 1;
end

% Initialize common variables
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

U = sqrt( nu(1)^2 + nu(2)^2 );          % Speed

%% Ship model parameters
L = 83;                     % Length (m)
B = 18;                     % Beam (m)
T = 5;                      % Draft (m)
rho = 1025;                 % Density of water (kg/m3)
Cb = 0.65;                  % Block coefficient: Cb = nabla / (L * B * T)
S = L * B + 2 * T * B;      % Wetted surface, box approximation

% Thrust: T = K_max * abs(n/n_max) * (n/n_max) = K_thr * abs(n) * n
K_max = [300e3 300e3 420e3 655e3]';        % Max propeller thrust (N)
n_max = [140 140 150 200]';                % Max propeller speed (rpm)
K_thr = diag(K_max./n_max.^2);             % Thruster coefficient matrix
l_x = [37, 35, -L/2, -L/2];                % Thruster x-coordinates
l_y = [0, 0, 7, -7];                       % Thruster y-coordinates

thrust_max = K_max(3)+K_max(4);     % max thrust in the surge direction (N)
U_max = 7.7;           % Max cruise speed (m/s) corresponding to max thrust

nabla = Cb * L * B * T;     % Volume displacement(m3)
m = rho * nabla;            % Mass (kg)
r_bg = [-4.5 0 -1.2]';      % Location of the CG with respect to the CO

Cw = 0.8;                   % Waterplane area coefficient: Cw = Awp/(L * B)
Awp = Cw * B * L;           % Waterplane area
KB = (1/3) * (5*T/2 - nabla/Awp);                         % Eq. (4.38)
k_munro_smith =  (6 * Cw^3) / ( (1+Cw) * (1+2*Cw));       % Eq. (4.37)
r_bb = [-4.5 0 T-KB]';      % Location of the CB with respect to the CO
BG = r_bb(3) - r_bg(3);     % Vertical distance between CG and CB

I_T = k_munro_smith * (B^3 * L) / 12;   % Transverse moment of inertia
I_L = 0.7 * (L^3 * B) / 12;             % Longitudinal moment of inertia
BM_T = I_T / nabla;
BM_L = I_L / nabla;
GM_T = BM_T - BG;                       % Should be larger than 0.5 m
GM_L = BM_L - BG;

% G matrix
LCF = -0.5;                   % x-distance from the CO to the center of Awp
r_bp = [0 0 0]';              % Compute G in the CO
G = Gmtrx(nabla,Awp,GM_T,GM_L,LCF,r_bp);

% Rigid-body mass matrix MRB
R44 = 0.35 * B;          % Radius of gyration in roll, see Eq.(4.77)-(4.78)
R55 = 0.25 * L;          % Radius of gyration in pitch
R66 = 0.25 * L;          % Radius of gyration in yaw
[MRB,CRB] = rbody(m,R44,R55,R66,nu(4:6),r_bg');     % MRB and CRB in the CG

% The added mass matrix MA is derived from the frequency-dependent potential 
% coefficients using a look-alike supply vessel in the MSS toolbox. The data
% is stored in the structure <vessel>
%   load supply.mat             % Check data by typing vessel
%   disp(vessel.main)
%   plotABC(vessel,'A')
%   w_0 = vessel.freqs;         % The minimum frequency w_0 is approximated
%   MA = vessel.A(:,:,1);       % as the zero-frequency used to compute MA 
MA = 1e10 * [   0.0001    0         0         0         0         0
                0    0.0003         0    0.0004         0   -0.0006
                0         0    0.0019         0    0.0233         0
                0    0.0004         0    0.0051         0   -0.0103
                0         0    0.0233         0    1.0816         0
                0   -0.0006         0   -0.0103         0    0.1156 ];

% Calibration using the values corresponding to the natural frequencies in
% heave, roll and pitch, found from visual inspection of the plots
MA(3,3) = 0.8e7;     % plotABC(vessel,'A',3,3,1)
MA(4,4) = 10.2e7;    % plotABC(vessel,'A',4,4,1)
MA(5,5) = 4.2e9;     % plotABC(vessel,'A',5,5,1)

CA  = m2c(MA, nu);   % Hydrodynamic added Coriolis and centripetal forces

% Mass matrix including hydrodynamic added mass
M = MRB + MA;
Minv = inv(M);

% D matrix
T1 = 100;               % Time constants for linear damping (s)
T2 = 100;
T6 = 1;
zeta4 = 0.15;           % Relative damping ratio in roll
zeta5 = 0.3;            % Relative damping ratio in pitch
D = Dmtrx([T1, T2, T6],[zeta4,zeta5],MRB,MA,G);

% Ocean current
v_c = [ Vc * cos(betaVc - eta(6))
        Vc * sin(betaVc - eta(6))
        0                          ];
nu_c = [v_c' zeros(1,3) ]';
nu_c_dot = [-Smtrx(nu(4:6)) * v_c
             zeros(3,1)           ];
nu_r = nu - nu_c;

% Add linear and quadratic drag in surge using the blending function sigma
[X,Xuu,Xu] = forceSurgeDamping(flag,nu_r(1),m,S,L,T1,rho,U_max,thrust_max);
tau_drag = [ X; zeros(5,1) ];

% Avoid double counting, linear and quadratic damping terms
D(1,1) = 0; % using: X = sigma * Xu * u_r + (1 - sigma) * Xuu * abs(u_r)*u_r

% Add crossflow drag
tau_crossflow = crossFlowDrag(L,B,T,nu_r);

% Thrust: T = K_max * u, where u = abs(n/n_max) * (n/n_max) is normalized
u_thr = abs(ui(1:4)) .* ui(1:4);               % Quadratic propeller speed
alpha = ui(5:6);                               % Azimuth angles
T_thr = thrConfig( {'T', 'T', alpha(1), alpha(2)}, l_x, l_y);
tau_3dof = T_thr * K_thr * u_thr;
tau_thr = [ tau_3dof(1) tau_3dof(2) 0 0 0 tau_3dof(3) ]';

% Kinematics
J = eulerang(eta(4),eta(5),eta(6));

% Equations of motion (Fossen 2021, Eqs. 6.111-6.116)
eta_dot = J * nu;
nu_dot = nu_c_dot + ...
    Minv * ( tau_thr +  tau_drag + tau_crossflow...
    - (CRB + CA + D) * nu_r - G * eta);

xdot = [nu_dot; eta_dot];    

%% Print vessel data
if nargin == 0 && nargout == 0

    % Natural frequencies
    w3 = sqrt( G(3,3) / M(3,3) );
    w4 = sqrt( G(4,4) / M(4,4) );
    w5 = sqrt( G(5,5) / M(5,5) );
    T3 = 2 * pi / w3;
    T4 = 2 * pi / w4;    
    T5 = 2 * pi / w5;

    fprintf('\n');
    fprintf('%s\n','-------------------------------------------------------------------------------------');
    fprintf('%s\n','OFFSHORE SUPPLY VESSEL MAIN CHARACTERISTICS');
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
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in roll:', zeta4);
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in pitch:', zeta5);
    fprintf('%-40s %8.2f m \n', 'Transverse metacentric height (GM_T):', GM_T);
    fprintf('%-40s %8.2f m \n', 'Longitudinal metacentric height (GM_L):', GM_L);
    fprintf('%-40s %8.2f s \n', 'Natural period in heave (T3):', T3);
    fprintf('%-40s %8.2f s \n', 'Natural period in roll (T4):', T4);
    fprintf('%-40s %8.2f s \n', 'Natural period in pitch (T5):', T5);   
    fprintf('%-40s %8.2f \n', 'Linear surge damping coefficient (Xu):', Xu); 
    fprintf('%-40s %8.2f \n', 'Quadratic drag coefficient (X|u|u):', Xuu); 
   
    D(1,1) = -Xu;
    matrices = {'Mass matrix: M = MRB + MA', M;...
        'Linear damping matrix: D', D; 'Restoring matrix: G', G};
    
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