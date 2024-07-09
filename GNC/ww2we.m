function Omega_e = ww2we(chi, U, Omega)
% Function to transform from wave frequency to encounter frequency
%
% Usage: Omega_e = ww2we(chi, U, Omega)
%
% Inputs:
%   chi   - Encounter angle [rad], 0 for following seas, pi for head seas
%   U     - Forward speed [m/sec]
%   Omega - Vector of wave frequency values [rad/sec]
%
% Outputs:
%   Omega_e - Vector of encounter frequency values [rad/sec]
%
% Reference: 
%   T. I. Fossen (2021) "Handbook of Marine Craft Hydrodynamics and Motion
%     Control, 2nd edtion, John Wiley & Sons Ltd., Chichester, UK.  
%
% Created by: Thor I. Fossen
% Date: 2024-07-09

g = 9.81; % Acceleration of gravity

% Calculate the encounter frequency 
Omega_e = Omega - (Omega.^2 * U * cos(chi) / g);

end