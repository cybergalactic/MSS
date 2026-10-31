function results = exRegularWaveEstimation
% exRegularWaveEstimation compares manual regular-wave excitation estimates
% with MSS-Capytaine force RAOs for a supply vessel and the Maritime Robotics
% Otter USV at 80 kg. Both craft are at zero forward speed in deep water.
%
% Manual estimates:
%   F3 = rho*g*Awp*A          (heave force in head seas)
%   F4 = rho*g*nabla*GM_T*k*A (roll moment in beam seas)
%   F5 = rho*g*nabla*GM_L*k*A (pitch moment in head seas)
% where A is wave amplitude, k = omega^2/g and omega = 2*pi/T.
%
% The MAT-files contain hydrostatic stiffness and force RAOs about the CG.
% Both are transformed to the center of flotation (CF) before comparison.
% Awp is recovered from the hydrostatic heave stiffness; the metacentric
% heights and displacement are read from vessel.main. No force RAOs are
% needed for the manual estimates. The six complex force RAOs are
% interpolated before shifting their moments from the CG to the CF.
%
% The Otter data use the approximate hull geometry included in MSS-Capytaine.
% This comparison illustrates the estimation methods; it is not a comparison
% with measured wave loads.
%
% Usage:
%   exRegularWaveEstimation           % Print the comparison tables
%   results = exRegularWaveEstimation; % Also return the computed values
%
% results contains the wave conditions, hydrostatic properties and the
% F_manual and F_RAO vectors, ordered as [F3; F4; F5], in N, Nm and Nm.
% The supply-vessel case uses A = 1 m and T = 20 s. This represents
% long-period swell: lambda = g*T^2/(2*pi) is approximately 625 m, about
% 7.5 times the vessel length of 82.8 m. The wave height is 2*A = 2 m.
% This long-wave case is chosen to illustrate the hydrostatic estimates.
% The Otter case uses A = 0.1 m and T = 3 s, giving a wavelength of 14.05 m.
% Change these settings below to examine other waves.
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK, Chapter 4.
%
% Author:    Thor I. Fossen
% Date:      2026-10-30
% Revisions:

% Load the two MSS-Capytaine craft and select regular-wave conditions.
load supply
craft(1).name = 'Supply vessel';
craft(1).vessel = vessel;
craft(1).A = 1.0;                 % Wave amplitude (m)
craft(1).T = 20.0;                % Wave period (s)

load otter_80kg
craft(2).name = 'Otter 80 kg';
craft(2).vessel = vessel;
craft(2).A = 0.1;                 % Wave amplitude (m)
craft(2).T = 3.0;                 % Wave period (s)

% Heave and pitch in head seas; roll in beam seas. MSS headings specify
% propagation direction: pi = head seas, pi/2 = beam seas.
dofs = [3 4 5];
beta_rel = [pi pi/2 pi];
loadNames = {'F3: heave force', 'F4: roll moment', 'F5: pitch moment'};
seaNames = {'Head', 'Beam', 'Head'};

for c = 1:numel(craft)
    vessel = craft(c).vessel;
    A = craft(c).A;
    T = craft(c).T;
    g = vessel.main.g;
    rho = vessel.main.rho;
    nabla = vessel.main.nabla;
    omega = 2*pi/T;
    k = omega^2/g;
    wavelength = 2*pi/k;

    speedIndex = find(abs(vessel.velocities) < 1e-12,1);
    if isempty(speedIndex)
        error('%s has no zero-speed hydrodynamic data.',craft(c).name);
    end
    w = vessel.forceRAO.w(:);
    if omega < min(w) || omega > max(w)
        error('Wave frequency for %s is outside the force-RAO grid.', ...
            craft(c).name);
    end

    % r_FG points from CF to CG. H maps CF velocity to CG velocity;
    % its transpose maps generalized loads from CG to CF.
    r_FG = vessel.main.CG(:) - vessel.main.CF(:);
    H = Hmtrx(r_FG);
    C_F = H' * vessel.C(:,:,1,speedIndex) * H;
    Awp = C_F(3,3)/(rho*g);
    GM_T = vessel.main.GM_T;
    GM_L = vessel.main.GM_L;

    % Manual estimates from the hydrostatic properties and wave slope k*A.
    F_manual = rho*g*A * [Awp; nabla*GM_T*k; nabla*GM_L*k];

    % Interpolate all six complex loads, then shift moments before taking
    % magnitudes. Heave couples into pitch when changing reference point.
    F_RAO = zeros(3,1);
    for n = 1:numel(dofs)
        headingError = abs(mod(vessel.headings-beta_rel(n)+pi,2*pi)-pi);
        [headingError,headingIndex] = min(headingError);
        if headingError > 1e-10
            error('The required wave heading is absent for %s.',craft(c).name);
        end
        H_G = zeros(6,1);
        for i = 1:6
            forceRAO = vessel.forceRAO.amp{i}(:,headingIndex,speedIndex) .* ...
                exp(1i*vessel.forceRAO.phase{i}(:,headingIndex,speedIndex));
            H_G(i) = interp1(w,forceRAO,omega,'linear');
        end
        H_F = H' * H_G;
        F_RAO(n) = A*abs(H_F(dofs(n)));
    end

    results(c).name = craft(c).name;
    results(c).waveAmplitude = A;
    results(c).wavePeriod = T;
    results(c).omega = omega;
    results(c).wavelength = wavelength;
    results(c).beta_rel = beta_rel;
    results(c).Awp = Awp;
    results(c).GM_T = GM_T;
    results(c).GM_L = GM_L;
    results(c).F_manual = F_manual;
    results(c).F_RAO = F_RAO;

    fprintf('\n%s: A = %.3g m, T = %.3g s, wavelength = %.3f m\n', ...
        craft(c).name,A,T,wavelength);
    fprintf('Reference point: CF; speed: 0 m/s; deep water\n');
    fprintf('Awp = %.6g m^2, GM_T = %.6g m, GM_L = %.6g m\n', ...
        Awp,GM_T,GM_L);
    fprintf('%-19s %-6s %-5s %15s %15s\n', ...
        'Load','Seas','Unit','Manual','Force RAO');
    units = {'N','Nm','Nm'};
    for n = 1:numel(dofs)
        fprintf('%-19s %-6s %-5s %15.6f %15.6f\n', ...
            loadNames{n},seaNames{n},units{n},F_manual(n),F_RAO(n));
    end
end

if nargout == 0
    clear results
end
end
