function Omega_e = encounter(chi, U, Omega)
% Function to transform from wave frequency to encounter frequency
%
% Usage: Omega_e = encounter(chi, U, Omega)
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
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Created by: Thor I. Fossen
% Date: 2024-07-09

constant = mssConstants();
g = constant.g;

% Calculate the encounter frequency 
Omega_e = Omega - (Omega.^2 * U * cos(chi) / g);

end
