function msi = HMmsi(a_z, w_e)
% msi = HMmsi(a_z,w_e) computes the Motion Sickness Incidence using the 
% method of O'Hanlon and McCauley (1974).
%
% Inputs:
%   a_z: Vector of vertical accelerations (m/s^2) in NED
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

% Vertical acceleration RMS-value 
a_z_rms = sqrt( mean(a_z.^2) );

% Constants
g = 9.81;  % Gravity (m/s²)
epsilon = 1e-6; % Small value to prevent log(0)

% Convert encounter frequency from rad/s to Hz
f_e = w_e / (2 * pi);

% Compute mu_MSI from O'Hanlon & McCauley model
mu_MSI = -0.819 + 2.32 * (log10(max(f_e, epsilon))).^2;

% Compute I factor
I = (-log10(a_z_rms / g) + mu_MSI) / (0.4 * sqrt(2));

% Compute MSI using the error function (erf)
msi = 100 * (0.5 * (1 - erf(I)));

end