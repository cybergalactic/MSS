function constant = mssConstants()
% mssConstants returns common physical constants used by MSS.
%
% Output:
%   constant.rho_water : Seawater density (kg/m^3)
%   constant.rho_air   : Air density (kg/m^3)
%   constant.g         : Standard gravity (m/s^2)
%
% Author: Thor I. Fossen
% Date: 2026-10-09

constant.rho_water = 1025;   % Seawater density (kg/m^3)
constant.rho_air = 1.225;    % Air density (kg/m^3)
constant.g = 9.81;           % Standard gravity (m/s^2)

end
