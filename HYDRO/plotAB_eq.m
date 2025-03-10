function plotAB_eq(vessel, mtrx, velno)
% plotAB_eq plots all elements A_ij or B_ij versus frequency and overlays 
% A_eq or B_eq as a red line. The equivalent added mass A_eq and damping B_eq 
% matrices are computed by integrating the frequency-dependent hydrodynamic 
% coefficients A(ω) and B(ω) using the wave spectrum S(ω) as a weighting 
% function, see computeManeuveringModel.m.
%
% Inputs:
%   vessel - Structure containing hydrodynamic matrices (A, B, A_eq, B_eq)
%   mtrx   - 'A' for added mass, 'B' for damping
%   velno  - (Optional) Velocity index. Default is 1 if omitted.
%
% Example:
%    plotAB_eq(vessel, 'A')      % Plots A_ij for velno = 1
%    plotAB_eq(vessel, 'B', 2)   % Plots B_ij for velno = 2
%
% Author:    Thor I. Fossen
% Date:      2025-03-10

% Select matrix type
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
    error('Invalid matrix type. Use ''A'' for added mass or ''B'' for damping.');
end

% Extract frequency data
freqs = vessel.freqs;
nfreq = length(freqs);

% Determine if speed is present in the dataset
sizeH = size(H);
hasVelocity = (length(sizeH) == 4); % True if H is 6×6×36×nvel

if hasVelocity
    nvel = sizeH(4); % Number of velocity cases
    if nargin < 3  % Default to velno = 1 if not provided
        velno = 1;
    end
    if velno > nvel || velno < 1
        error('Invalid velocity index: velno must be in the range 1 to %d', nvel);
    end
    H = H(:,:,:,velno);  % Extract for selected velno
else
    velno = 1; % Default to 1 if no speed dependency
end

% Plot longitudinal elements (1,1), (1,3), (3,1), (3,3), etc.
figure(figno);
clf;
k = 1;
for i = 1:2:5
    for j = 1:2:5
        subplot(3,3,k);
        Hplot = reshape(H(i,j,:), nfreq, 1);
        plot(freqs, Hplot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(freqs, repmat(H_eq(i,j), size(freqs)), 'r-', 'LineWidth', 2); % Red line
        grid on;
        
        Hw = sprintf('%s_{%d%d}', mtrx, i, j);
        xlabel('Frequency (rad/s)');
        title(sprintf('%s - %s_{%d%d} (Vel %d)', titleStr, mtrx, i, j, velno));
        k = k + 1;
    end
end

% Plot lateral elements (2,2), (2,4), (4,2), (4,4), etc.
figure(figno + 1);
clf;
k = 1;
for i = 2:2:6
    for j = 2:2:6
        subplot(3,3,k);
        Hplot = reshape(H(i,j,:), nfreq, 1);
        plot(freqs, Hplot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(freqs, repmat(H_eq(i,j), size(freqs)), 'r-', 'LineWidth', 2); % Red line
        grid on;
        
        Hw = sprintf('%s_{%d%d}', mtrx, i, j);
        xlabel('Frequency (rad/s)');
        title(sprintf('%s - %s_{%d%d} (Vel %d)', titleStr, mtrx, i, j, velno));
        k = k + 1;
    end
end

end