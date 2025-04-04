function plotAB_eq(vessel, mtrx, velNo)
% plotAB_eq plots all elements A_ij or B_ij versus frequency ω = vessel.freqs 
% and overlays A_eq or B_eq for each wave spectrum peak frequency ω_p = vessel.omega_p. 
% The equivalent added mass A_eq and damping B_eq matrices are computed by 
% integrating the frequency-dependent hydrodynamic coefficients A(ω) and B(ω) 
% using the wave spectrum S(ω, ω_p) as a weighting function, 
% see computeManeuveringModel.m.

% Inputs:
%   vessel - Structure with fields A, B, A_eq, B_eq, freqs, omega_p, velocities
%   mtrx   - 'A' or 'B' to select added mass or damping matrices
%   velNo  - Velocity index (1-based). Default is 1.
%
% Author: Thor I. Fossen
% Date:   2025-03-10

% Handle matrix selection
if strcmp(mtrx, 'A')
    H = vessel.A;
    H_eq = vessel.A_eq;
    figno = 100;
    titleStr = 'Added Mass';
elseif strcmp(mtrx, 'B')
    H = vessel.B;
    H_eq = vessel.B_eq;
    figno = 200;
    titleStr = 'Damping';
else
    error('Invalid matrix type. Use ''A'' or ''B''.');
end

% Handle velocity index
if nargin < 3
    velNo = 1;
end

% Check velocity bounds
if isfield(vessel, 'velocities')
    nvel = length(vessel.velocities);
else
    nvel = 1;
end

if velNo < 1 || velNo > nvel
    error('Invalid velocity index: must be between 1 and %d', nvel);
end

% Frequency and omega_p vectors
freqs = vessel.freqs(:);
omega_p = vessel.omega_p(:);

% Extract 3D slice for current velocity
H_w = H(:,:,:,velNo);
H_eq = H_eq(:,:,:,velNo);  % [6×6×nOmega]

% ---- LONGITUDINAL ELEMENTS ----
figure(figno);
clf;
k = 1;
for i = 1:2:5
    for j = 1:2:5
        subplot(3,3,k);
        
        % Frequency-dependent matrix element and interpolation onto omega_p
        Hij = squeeze(H_w(i,j,:));

        % Plot interpolated H_w at omega_p
        plot(freqs, Hij, 'b-o', 'LineWidth', 1.5); hold on;

        % Plot equivalent H_eq values
        Hij_eq = squeeze(H_eq(i,j,:));
        plot(omega_p, Hij_eq, 'r-s', 'LineWidth', 1.5);

        xlabel('\omega (rad/s)');
        ylabel(sprintf('%s_{%d%d}', mtrx, i, j));
        title(sprintf('%s - %s_{%d%d} (Vel %d)', titleStr, mtrx, i, j, velNo));
        grid on;
        legend('Raw data', 'Equivalent', 'Location', 'best');
        k = k + 1;
    end
end

% ---- LATERAL ELEMENTS ----
figure(figno + 1);
clf;
k = 1;
for i = 2:2:6
    for j = 2:2:6
        subplot(3,3,k);
        
        Hij= squeeze(H_w(i,j,:));

        plot(freqs, Hij, 'b-o', 'LineWidth', 1.5); hold on;

        Hij_eq = squeeze(H_eq(i,j,:));
        plot(omega_p, Hij_eq, 'r-s', 'LineWidth', 1.5);

        xlabel('\omega (rad/s)');
        ylabel(sprintf('%s_{%d%d}', mtrx, i, j));
        title(sprintf('%s - %s_{%d%d} (Vel %d)', titleStr, mtrx, i, j, velNo));
        grid on;
        legend('Raw data', 'Equivalent', 'Location', 'best');
        k = k + 1;
    end
end

end