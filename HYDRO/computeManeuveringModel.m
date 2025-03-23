function vessel = computeManeuveringModel(vessel, velno, omega_p, plotFlag)
% computeManeuveringModel is compatible with MATLAB and GNU Octave (www.octave.org).
% Computes the equivalent added mass A_eq and damping B_eq matrices by
% integrating the frequency-dependent hydrodynamic coefficients A(ω) and B(ω) 
% using the wave spectrum S(ω) as a weighting function.
%
% The equivalent matrices are computed as:
%
%   A_eq(i,j) = ∫ A(i,j,ω) S(ω) dω  /  ∫ S(ω) dω
%   B_eq(i,j) = ∫ B(i,j,ω) S(ω) dω  /  ∫ S(ω) dω
%
% where:
%   - A(i,j,ω) and B(i,j,ω) are the added mass and damping coefficients 
%     at frequency ω for each matrix element (i,j).
%   - S(ω) is the wave energy spectrum, computed using waveSpectrum().
%   - The integrals are evaluated numerically using trapezoidal integration.
%
% This ensures that the kinetic energy and power dissipation properties 
% of the frequency-dependent system are preserved in the equivalent 
% constant matrices.
%
% Inputs:
%   vessel       - Structure containing vessel hydrodynamic data (A, B, freqs)
%   velno        - Velocity index (if speeds exist in vessel). If no speeds, velno = 1
%   omega_p      - Wave peak frequency for scaling of the model
%   plotFlag     - Set to 1 to plot the A(ω) and B(ω) matrix elememts, 0 for no plot
%
% Outputs:
%   vessel.A_eq - Equivalent added mass matrix
%   vessel.B_eq - Equivalent damping matrix
%
% Example call: 
%   load supply % other vessels: s175, tanker, fpso, semisub
%   w_p = 1.2;
%   vessel = computeManeuveringModel(vessel, 1, w_p, 1)
%   disp(vessel.A_eq)
%   disp(vessel.B_eq)
%
% Author: Thor I. Fossen
% Date: 2025-03-10

% Determine if speed is present in the matrices
sizeA = size(vessel.A);    % Get size of A matrix
hasVelocity = (length(sizeA) == 4);  % True if A is 6×6×36×nvel

if hasVelocity
    nvel = sizeA(4);  % Number of velocity cases
    if velno > nvel || velno < 1
        error('Invalid velocity index: velno must be in the range 1 to %d', nvel);
    end
    A_w = vessel.A(:,:,:,velno);  % Extract for selected velno
    B_w = vessel.B(:,:,:,velno); 
else
    disp('Only one velocity in the data set');
    velno = 1;  % Default to 1 if no speed dependency
    A_w = vessel.A;  % No speed dependency
    B_w = vessel.B;
end

% Define finer frequency grid for interpolation
freqs_fine = linspace(min(vessel.freqs), max(vessel.freqs), 100)'; 

% Compute the PM wave spectrum at given frequencies
alpha = 8.1e-3 * (9.81)^2;
beta = 0.74;
S_w = alpha ./ freqs_fine.^5 .* exp(-beta * (omega_p ./ freqs_fine).^4);

% Compute normalization factor (integral of S)
denom = trapz(freqs_fine, S_w);

% Initialize equivalent matrices
A_eq = zeros(6,6);
B_eq = zeros(6,6);

% Integrate each element of A and B separately
for i = 1:6
    for j = 1:6
        % Extract the (i,j) element for all frequencies
        A_ij_w = squeeze(A_w(i,j,:)); % Ensure it's a column vector
        B_ij_w = squeeze(B_w(i,j,:));

        % Interpolate A_w and B_w onto the finer grid
        A_interp = interp1(vessel.freqs, A_ij_w, freqs_fine, 'pchip'); 
        B_interp = interp1(vessel.freqs, B_ij_w, freqs_fine, 'pchip');

        % Numerically integrate element-wise using trapezoidal rule
        A_eq(i,j) = trapz(freqs_fine, A_interp .* S_w) / denom;
        B_eq(i,j) = trapz(freqs_fine, B_interp .* S_w) / denom;
    end
end

% Store results in vessel structure
vessel.A_eq = A_eq;
vessel.B_eq = B_eq;

% Plots
if plotFlag == 1
    plotAB_eq(vessel, 'A', velno);
    plotAB_eq(vessel, 'B', velno);
end

end