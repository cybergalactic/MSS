function vessel = computeManeuveringModel(vessel, omega_p, plotFlag)
% computeManeuveringModel is compatible with MATLAB and GNU Octave (www.octave.org).
% Computes the equivalent added mass vessel.A_eq(:,:,k,n) and damping
% vessel.B_eq(:,:,k,n) matrices for each omega_p(k) and velocity(n) by integrating
% the frequency-dependent hydrodynamic coefficients A(ω) and B(ω) using the wave
% spectrum S(ω,ω_p) as a weighting function.
%
% The equivalent matrices are computed as:
%
%   A_eq(i,j,ω_p,velocity) = ∫ A(i,j,ω,velocity) S(ω,ω_p(k)) dω / ∫ S(ω,ω_p(k)) dω
%   B_eq(i,j,ω_p,velocity) = ∫ B(i,j,ω,velocity) S(ω,ω_p(k)) dω / ∫ S(ω,ω_p(k)) dω
%
% where:
%   - A(i,j,ω,velocity) and B(i,j,ω,velocity) are the added mass and damping
%     coefficients at frequency and velocity for each matrix element (i,j).
%   - S(ω,ω_p) is the wave energy spectrum.
%   - The integrals are evaluated numerically using trapezoidal integration.
%
% This ensures that the kinetic energy and power dissipation properties
% of the frequency-dependent system are preserved in the equivalent
% constant matrices.
%
% Inputs:
%   vessel       - Structure containing vessel hydrodynamic data (A, B, freqs)
%   omega_p      - Vector of wave peak frequencies, e.g. omega_p = linspace(0.1, 3.0, 10)
%   plotFlag     - Set to 1 to plot the A(ω) and B(ω) matrix elememts, 0 for no plot
%
% Outputs:
%   vessel.omega_p - Vector of wave spectrum peak frequencies
%   vessel.A_eq    - Equivalent added mass matrix
%   vessel.B_eq    - Equivalent damping matrix
%
% Example call:
%   load supply % other vessels: s175, tanker, fpso, semisub
%   omega_p = linspace(0.1, 3.0, 20)
%   vessel = computeManeuveringModel(vessel, 1, omega_p, 1)
%   disp(vessel.A_eq)
%   disp(vessel.B_eq)
%
% Author: Thor I. Fossen
% Date: 2025-03-10

% Check number of velocity cases
if isfield(vessel, 'velocities') && ~isempty(vessel.velocities)
    nvel = length(vessel.velocities);
else
    nvel = 1;
end

A_all = vessel.A; % Added mass
B_all = vessel.B + vessel.Bv; % Total Damping, sum of potential and viscous damping 

omega_min = min(vessel.freqs);
omega_max = max(vessel.freqs);

% Avoid omega = 0 to prevent numerical issues in spectrum normalization
if omega_min == 0
    omega_min = 1e-6;
end

% Define finer frequency grid for interpolation
freqs_fine = linspace(omega_min, omega_max, 100)';
nOmega = length(omega_p);

% PM wave spectrum parameters
alpha = 8.1e-3 * (9.81)^2;
beta = 0.74;

% Initialize equivalent matrices of dimension [6, 6, nOmega, nvel]
Aeq_all = zeros(6, 6, nOmega, nvel);
Beq_all = zeros(6, 6, nOmega, nvel);

% Loop over all velocities
for velNo = 1:nvel
    A_w = A_all(:,:,:,velNo);
    B_w = B_all(:,:,:,velNo);

    % Loop over all omega_p
    for k = 1:nOmega
        op = omega_p(k);
        S_w = alpha ./ freqs_fine.^5 .* exp(-beta * (op ./ freqs_fine).^4);
        denom = trapz(freqs_fine, S_w);

        % Loop over DOFs
        for i = 1:6
            for j = 1:6
                A_ij_w = squeeze(A_w(i,j,:));
                B_ij_w = squeeze(B_w(i,j,:));

                A_interp = interp1(vessel.freqs, A_ij_w, freqs_fine, 'pchip');
                B_interp = interp1(vessel.freqs, B_ij_w, freqs_fine, 'pchip');

                Aeq_all(i,j,k,velNo) = trapz(freqs_fine, A_interp .* S_w) / denom;
                Beq_all(i,j,k,velNo) = trapz(freqs_fine, B_interp .* S_w) / denom;
            end
        end
    end
end

% Store in vessel
vessel.omega_p = omega_p;
vessel.A_eq = Aeq_all;
vessel.B_eq = Beq_all;

% Optional plotting for velcoity #1
if plotFlag == 1
    plotAB_eq(vessel, 'A', 1);
    plotAB_eq(vessel, 'B', 1);
end

end