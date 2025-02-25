function msi = HMmsi(a_z, w_e)
% msi = HMmsi(a_z,w_e) computes the Motion Sickness Incidence using the 
% method of O'Hanlon and McCauley (1974).
%
% Inputs:
%   a_z: RMS of vertical acceleration (m/s^2) expressed in NED
%   w_e: Vector of encounter frequencies (rad/s)
%
% Outputs:
%   msi: The percentage of persons that become seasick during a 2-hour voyage
%
% References:
%   - A. R. J. M. Lloyd (1989). Seakeeping Behaviour in Rough Water.
%   - E. V. Lewis (1989). Principles of Naval Architecture. Vol III.
%   - J. F. O'Hanlon and M. E. McCauley (1974). Aerospace Medicine.
%
% Author:    J. M. de la Cruz
% Date:      2000-03-28
% Revisions: 
%   2024-08-07 : Minor changes of notation/documentation. 
%   2025-02-25 : Replaced abs(a_z) with the rms value. Converted encounter 
%                frequencies from rad/s to Hz

% Ensure inputs are column vectors
a_z = a_z(:);
w_e = w_e(:);

g = 9.81;  % Acceleration due to gravity (m/s^2)
epsilon = 1e-10; % Avoid log(0) by adding a small epsilon

% Convert encounter frequency from rad/s to Hz
f_e = w_e / (2 * pi);  

% Compute mu_MSI based on O'Hanlon and McCauley's model
mu_MSI = -0.819 + 2.32 * (log10(f_e + epsilon)).^2; 

% Compute I factor for each frequency
I = ( -log10(a_z / g) + mu_MSI ) / 0.4;

% Divide by sqrt(2) to adjust for MATLAB's erf function
I = I / sqrt(2);

% Compute MSI for each frequency
msi = zeros(length(w_e), 1);
for i = 1:length(w_e)
    if I(i) >= 0
        msi(i) = 0.5 - 0.5 * erf(I(i));
    else
        msi(i) = 0.5 + 0.5 * erf(-I(i));  
    end
end

% Convert to percentage
msi = 100 * msi;