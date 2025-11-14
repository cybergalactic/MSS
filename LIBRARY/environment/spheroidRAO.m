function vessel = spheroidRAO(vessel,a,b,zn,verbose)
% spheroidRAO computes 6-DOF *first-order Froude–Krylov (FK)* wave excitation 
% forces and moments for a submerged prolate spheroid.
%
% The output is stored in vessel.forceRAO.* using the same fields as the MSS 
% toolbox (Re, Im, amp, phase), enabling direct use for time-domain wave 
% reconstruction.
%
% This model includes ONLY the incident-wave pressure field (FK excitation) 
% and does not include diffraction or radiation effects.
%
% INPUTS:
%   vessel : existing vessel structure or [] to create a new structure
%   a      : spheroid semi-major axis (a > b), half-length     [m]
%   b      : spheroid semi-minor axis (radius)                 [m]
%   zn     : submergence depth of the spheroid center (> 0)    [m]
%   verbose: true → plot FK curves for 0–180° headings (optional)
%
% OUTPUT (stored in vessel.forceRAO):
%   vessel.forceRAO.w        : wave frequencies ω [rad/s]
%   vessel.forceRAO.Re{dof}  : real(FK)   [N or Nm]
%   vessel.forceRAO.Im{dof}  : imag(FK)=0
%   vessel.forceRAO.amp{dof} : |FK(ω,β)|
%   vessel.forceRAO.phase{dof}: FK phase = 0
%   vessel.headings : wave headings from 0° to 350° in 10° steps, stored in
%                     stored in radians (36 unique directions).
%
% FK is 6 × Nω × Nβ:
%   DOFs: 1=surge, 2=sway, 3=heave, 4=roll, 5=pitch, 6=yaw
%
% ------------------------------------------------------------------------------
% THEORY — FIRST-ORDER FROUDE–KRYLOV EXCITATION
% ------------------------------------------------------------------------------
% The incident-wave dynamic pressure is:
%
%   p_FK = ρ g A exp(k z) cos(k ξ − ω t)
%
% where:
%   k  = ω²/g   (deep water)
%   z  = vertical coordinate (negative below surface)
%   ξ  = x cosβ + y sinβ (wave propagation direction)
%
% For a submerged body at depth zn:
%       exp(k z) = exp(−k zn)
%
% For an axisymmetric prolate spheroid, the generalized FK loads become:
%
%   Heave force: F3 = ρ g π b² exp(−k zn)
%   Surge:      F1 =  F3 cos(β)
%   Sway:       F2 =  F3 sin(β)
%
%   Pitch moment:  M5 = F3 * r_pitch * k
%       r_pitch = (a² − b²)/(5a)
%
%   Roll moment:   M4 = M5 sin(β)
%   Pitch moment:  M5 = M5 cos(β)
%
%   Yaw moment:    M6 = 0
%
% These form FK(ω, β) which is stored as the RAO table.
% ------------------------------------------------------------------------------
% EXAMPLES:
%   vessel = spheroidRAO(vessel,a,b,zn);
%   vessel = spheroidRAO([],2,1,5,true);
%
% Reference:
%   Fossen, T. I. (2021). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 2nd edtion. John Wiley & Sons Ltd., Chichester, UK.
%
% Author:       T.I. Fossen
% Data:         2025-11-12
% Revisions:

if nargin < 4
    verbose = false;
end

% Physical constants
g   = 9.81;
rho = 1025;

if zn <= 0
    error('zn must be positive (submerged)');
end

% ------------------------------------------------------------------------------
% Frequency vector (deep-water approximation)
% ------------------------------------------------------------------------------
omega = linspace(0.05,4,100)';     % rad/s
k = omega.^2 / g;                  % deep water
decay = exp(-k * zn);              % FK attenuation

% ------------------------------------------------------------------------------
% FK excitation formulas for a submerged spheroid
% ------------------------------------------------------------------------------
F3 = rho * g * pi * b^2 .* decay;      % Heave
F1 = F3;                                % Surge projection handled in β loop

r_pitch = (a^2 - b^2) / (5*a);
F5 = r_pitch .* k .* F3;                % Pitch moment

% ------------------------------------------------------------------------------
% Wave headings 0–180°, then mirror to 360°
% ------------------------------------------------------------------------------
beta_half = 0:10:180;
beta_rad  = deg2rad(beta_half);
Nbeta_half = length(beta_rad);
Nomega = length(omega);

% Preallocate FK tables: 6 × Nω × Nβ
FK_half = zeros(6,Nomega,Nbeta_half);

for ib = 1:Nbeta_half
    br = beta_rad(ib);

    FK_half(:, :, ib) = [
        (F1 .* cos(br))'                 % Surge
        (F1 .* sin(br))'                 % Sway
        F3'                              % Heave
        (F5 .* sin(br))'                 % Roll moment
        (F5 .* cos(br))'                 % Pitch moment
        zeros(1,Nomega)                  % Yaw moment
    ];
end


% ------------------------------------------------------------------------------
% Mirror FK from 0–180° to 180–360° using BODY-frame symmetry
% ------------------------------------------------------------------------------
FK_full = zeros(6, Nomega, 36);

% Copy 0–180° (19 headings)
FK_full(:,:,1:19) = FK_half;

% Mirror 10°–170° → 190°–350°
src = 18:-1:2;   % indices 2..18 reversed

% BODY-axis parity for FK forces/moments
% DOFs: [surge sway heave roll pitch yaw]
fk_sign = [1  -1  1  -1  1  -1];

for dof = 1:6
    FK_full(dof,:,20:36) = fk_sign(dof) * FK_half(dof,:,src);
end

for dof = 1:6
    vessel.forceRAO.Re{dof} = squeeze(FK_full(dof,:,:));
    vessel.forceRAO.Im{dof} = zeros(size(vessel.forceRAO.Re{dof}));
    vessel.forceRAO.amp{dof} = abs(vessel.forceRAO.Re{dof});
    vessel.forceRAO.phase{dof} = zeros(size(vessel.forceRAO.Re{dof}));
end

% ------------------------------------------------------------------------------
% Populate MSS-compatible RAO structure
% ------------------------------------------------------------------------------
vessel.forceRAO.w = omega;
vessel.headings = deg2rad(0:10:350)';    % 36 headings (stored in radians)

for dof = 1:6
    H = squeeze(FK_full(dof,:,:));   % [Nω × Nβ]
    vessel.forceRAO.Re{dof}    = H;
    vessel.forceRAO.Im{dof}    = zeros(size(H));
    vessel.forceRAO.amp{dof}   = abs(H);
    vessel.forceRAO.phase{dof} = zeros(size(H));
end

% ------------------------------------------------------------------------------
% OPTIONAL PLOTTING
% ------------------------------------------------------------------------------
if verbose
    labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
    figure; clf;
    for dof = 1:6
        subplot(3,2,dof); hold on; grid on;
        for ib = 1:length(beta_half)
            plot(omega, squeeze(FK_half(dof,:,ib)), 'LineWidth', 1.2);
        end
        title(labels{dof});
        xlabel('\omega [rad/s]');
        ylabel('FK amplitude');
    end
    sgtitle('Froude–Krylov Excitation Forces and Moments (0–180°)');
end

end
