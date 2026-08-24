function vessel = spheroidRAO(vessel,a,b,zn,verbose,includeDiffraction)
% spheroidRAO computes 6-DOF first-order wave excitation loads for a fully
% submerged prolate spheroid.
%
% The output is stored in vessel.forceRAO.* using the same fields as the MSS 
% toolbox (Re, Im, amp, phase), enabling direct use for time-domain wave 
% reconstruction.
%
% The Froude-Krylov (FK) loads are computed from the exact incident-pressure
% integral. Diffraction can be included using a low-frequency (LF) added-mass
% approximation. Radiation effects are not included.
%
% INPUTS:
%   vessel : existing vessel structure or [] to create a new structure
%   a      : spheroid semi-major axis (a > b), half-length     [m]
%   b      : spheroid semi-minor axis (radius)                 [m]
%   zn     : depth of the spheroid center (zn > b)             [m]
%   verbose: true -> plot excitation curves for 0-180 deg headings (optional)
%   includeDiffraction: true -> add LF diffraction forces (default: true)
%
% OUTPUT (stored in vessel.forceRAO):
%   vessel.forceRAO.w          : wave frequencies omega [rad/s]
%   vessel.forceRAO.Re{dof}    : real(excitation)  [N/m or Nm/m]
%   vessel.forceRAO.Im{dof}    : imag(excitation)  [N/m or Nm/m]
%   vessel.forceRAO.amp{dof}   : abs(excitation)   [N/m or Nm/m]
%   vessel.forceRAO.phase{dof} : angle(excitation) [rad]
%   vessel.headings            : 0-350 deg in 10 deg steps [rad]
%   vessel.forceRAO_FK          : FK contribution
%   vessel.forceRAO_diffraction : LF diffraction contribution
%
% DOFs: 1=surge, 2=sway, 3=heave, 4=roll, 5=pitch, 6=yaw.
%
% ------------------------------------------------------------------------------
% THEORY - FIRST-ORDER FROUDE-KRYLOV EXCITATION
% ------------------------------------------------------------------------------
% The incident-wave dynamic pressure per unit wave amplitude is
%
%   p = real(p_hat*exp(i*omega*t)),
%   p_hat = rho*g*exp(k*z)*exp(-i*k*xi),
%
% where k = omega^2/g, z is positive upward, and
%
%   xi = x*cos(beta) + y*sin(beta).
%
% The complex FK load on a closed body follows from the divergence theorem:
%
%   F_hat = -integral_S(p_hat*n)dS = -integral_V(grad(p_hat))dV.
%
% For x^2/a^2 + y^2/b^2 + z^2/b^2 <= 1, define
%
%   e      = sqrt(a^2-b^2),        lambda = k*e*cos(beta),
%   Phi    = 3*(sin(lambda)-lambda*cos(lambda))/lambda^3,
%   Psi    = 3*((3-lambda^2)*sin(lambda)-3*lambda*cos(lambda))/lambda^5,
%   V      = 4*pi*a*b^2/3,
%   P      = rho*g*k*V*Phi*exp(-k*zn),
%   Q      = rho*g*k^2*V*e^2*Psi*exp(-k*zn).
%
% The removable limits are Phi(0)=1 and Psi(0)=1/5. In the MSS body-fixed
% NED convention, the generalized FK load used by waveForceRAO is
%
%   [X Y Z K M N]' = [i*P*cos(beta), i*P*sin(beta), P, ...
%                      0, i*Q*cos(beta), Q*cos(beta)*sin(beta)]'.
%
% In the LF limit, the incident-wave acceleration is nearly uniform over the
% body. The diffraction force is then approximated by
%
%   F_diff = M_A * a_wave,
%
% where M_A is the unbounded-fluid translational added-mass matrix of the
% prolate spheroid. This correction applies to surge, sway, and heave. The
% rotational diffraction loads are zero at this order; finite-wavelength
% diffraction moments require a full diffraction solution. With
% P0 = rho*g*k*V*exp(-k*zn), the implemented diffraction RAO is
%
%   [X_D Y_D Z_D K_D M_D N_D]' = ...
%       [i*CAx*P0*cos(beta), i*CAt*P0*sin(beta), CAt*P0, 0, 0, 0]'.
%
% The approximation is intended for ka << 1 and sufficient submergence that
% free-surface corrections to the unbounded-fluid added mass are small.
% ------------------------------------------------------------------------------
% EXAMPLES:
%   vessel = spheroidRAO(vessel,a,b,zn);
%   vessel = spheroidRAO([],2,1,5,true);           
%   vessel = spheroidRAO([],2,1,5,true,false);     % FK only
%
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%       Control, 2nd edition. John Wiley & Sons Ltd., Chichester, UK.
%   Imlay, F. H. (1961). The Complete Expressions for Added Mass of a Rigid
%       Body Moving in an Ideal Fluid. DTMB Report 1528.
%
% Author:       T.I. Fossen
% Date:         2025-11-12
% Revisions:    2026-08-24 Corrected closed-surface FK loads and RAO phases.
%                          Added optional LF added-mass diffraction forces.

if nargin < 5 || isempty(verbose)
    verbose = false;
end
if nargin < 6 || isempty(includeDiffraction)
    includeDiffraction = true;
end

if isempty(vessel)
    vessel = struct();
elseif ~isstruct(vessel)
    error('vessel must be a structure or empty');
end

if ~isPositiveFiniteScalar(a) || ~isPositiveFiniteScalar(b) || a <= b
    error('a and b must be finite positive scalars satisfying a > b');
end
if ~isPositiveFiniteScalar(zn) || zn <= b
    error('zn must be a finite scalar satisfying zn > b (fully submerged)');
end
if ~isscalar(verbose)
    error('verbose must be a scalar');
end
if ~isscalar(includeDiffraction)
    error('includeDiffraction must be a scalar');
end

% Physical constants
if isfield(vessel,'main') && isfield(vessel.main,'g') && ...
        isPositiveFiniteScalar(vessel.main.g)
    g = vessel.main.g;
else
    g = 9.81;
end
rho = 1025;

% Deep-water wave frequencies and headings
omega = linspace(0.05,4,100)';
k = omega.^2 / g;
decay = exp(-k * zn);
beta = deg2rad(0:10:350);

% Spheroid geometry
volume = 4 * pi * a * b^2 / 3;
e2 = a^2 - b^2;
e = sqrt(e2);

% Unbounded-fluid translational added masses. The ellipsoidal potential
% coefficients satisfy alphaX + 2*alphaT = 2, and the added-mass ratios are
% C_A = alpha/(2-alpha). Here T denotes either transverse direction.
[CAx,CAt] = prolateAddedMassCoefficients(a,b);
addedMass = rho * volume * diag([CAx CAt CAt]);

Nomega = length(omega);
Nbeta = length(beta);
FK = complex(zeros(6,Nomega,Nbeta));
diffraction = complex(zeros(6,Nomega,Nbeta));

for ib = 1:Nbeta
    cb = cos(beta(ib));
    sb = sin(beta(ib));
    lambda = k * e * cb;

    Phi = spheroidForceFactor(lambda);
    Psi = spheroidMomentFactor(lambda);
    P = rho * g * k .* decay * volume .* Phi;
    Q = rho * g * k.^2 .* decay * volume * e2 .* Psi;
    P0 = rho * g * k .* decay * volume;

    FK(:,:,ib) = [
        (1i * P * cb).'             % Surge
        (1i * P * sb).'             % Sway
        P.'                          % Heave
        zeros(1,Nomega)              % Roll moment
        (1i * Q * cb).'             % Pitch moment
        (Q * cb * sb).'              % Yaw moment
    ];

    diffraction(:,:,ib) = [
        (1i * CAx * P0 * cb).'       % Surge
        (1i * CAt * P0 * sb).'       % Sway
        (CAt * P0).'                 % Heave
        zeros(1,Nomega)              % Roll moment
        zeros(1,Nomega)              % Pitch moment
        zeros(1,Nomega)              % Yaw moment
    ];
end

if includeDiffraction
    excitation = FK + diffraction;
    diffractionModel = 'low-frequency added-mass approximation';
else
    excitation = FK;
    diffractionModel = 'none';
end

% Populate MSS-compatible RAO structure
vessel.forceRAO.w = omega;
vessel.headings = beta';
vessel.forceRAO.diffractionModel = diffractionModel;
vessel.forceRAO.addedMass = addedMass;
vessel.forceRAO_FK = struct();
vessel.forceRAO_diffraction = struct();

for dof = 1:6
    H = squeeze(excitation(dof,:,:)); % [Nomega x Nbeta]
    H_FK = squeeze(FK(dof,:,:));
    H_diff = squeeze(diffraction(dof,:,:));

    vessel.forceRAO = storeComplexRAO(vessel.forceRAO,H,dof);
    vessel.forceRAO_FK = storeComplexRAO(vessel.forceRAO_FK,H_FK,dof);
    vessel.forceRAO_diffraction = storeComplexRAO( ...
        vessel.forceRAO_diffraction,H_diff,dof);
end
vessel.forceRAO_FK.w = omega;
vessel.forceRAO_diffraction.w = omega;

% Optional plotting
if verbose
    labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
    figure; clf;
    for dof = 1:6
        subplot(3,2,dof); hold on; grid on;
        for ib = 1:19
            plot(omega, abs(squeeze(excitation(dof,:,ib))), 'LineWidth', 1.2);
        end
        title(labels{dof});
        xlabel('\omega [rad/s]');
        ylabel('Excitation amplitude');
    end
    if ~isoctave
        if includeDiffraction
            sgtitle('FK + LF Diffraction Excitation Loads (0-180 deg)');
        else
            sgtitle('Froude-Krylov Excitation Loads (0-180 deg)');
        end
    end
end

end

function tf = isPositiveFiniteScalar(value)
tf = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
end

function RAO = storeComplexRAO(RAO,H,dof)
RAO.Re{dof} = real(H);
RAO.Im{dof} = imag(H);
RAO.amp{dof} = abs(H);
RAO.phase{dof} = angle(H);
end

function [CAx,CAt] = prolateAddedMassCoefficients(a,b)
% Added-mass ratios relative to displaced mass for a prolate spheroid.
eccentricity = sqrt(1 - (b/a)^2);
if eccentricity < 0.01
    e2 = eccentricity^2;
    alphaX = 2/3 - 4*e2/15 - 4*e2^2/35 - 4*e2^3/63;
else
    alphaX = (1 - eccentricity^2) / eccentricity^3 * ...
        (log((1 + eccentricity) / (1 - eccentricity)) - 2*eccentricity);
end
alphaT = 1 - alphaX / 2;
CAx = alphaX / (2 - alphaX);
CAt = alphaT / (2 - alphaT);
end

function Phi = spheroidForceFactor(lambda)
% Stable evaluation of 3*(sin(lambda)-lambda*cos(lambda))/lambda^3.
Phi = zeros(size(lambda));
small = abs(lambda) < 0.05;
x = lambda(small);
Phi(small) = 1 - x.^2/10 + x.^4/280 - x.^6/15120;
x = lambda(~small);
Phi(~small) = 3 * (sin(x) - x .* cos(x)) ./ x.^3;
end

function Psi = spheroidMomentFactor(lambda)
% Stable evaluation of the moment factor, whose limit at zero is 1/5.
Psi = zeros(size(lambda));
small = abs(lambda) < 0.05;
x = lambda(small);
Psi(small) = 1/5 - x.^2/70 + x.^4/2520 - x.^6/166320;
x = lambda(~small);
Psi(~small) = 3 * ((3 - x.^2) .* sin(x) - 3*x .* cos(x)) ./ x.^5;
end
