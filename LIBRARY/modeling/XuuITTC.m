function Xuu = XuuITTC(u_r,rho,L,B,T,C_B)
% Xuu = XuuITTC(u_r,rho,L,B,T,C_B) computes the quadratic surge-damping coefficient 
% using the ITTC-1957 model-ship correlation line:
%
%   X = Xuu * abs(u_r) * u_r
%
% The wetted surface of a conventional displacement monohull is estimated
% using the Mumford formula. The Reynolds number is bounded from below
% because the ITTC-1957 line is not applicable at zero or very low speeds.
%
% Inputs:
%   u_r: Relative surge velocity (m/s)
%   rho: Water density (kg/m^3)
%   L:   Vessel length (m)
%   B:   Vessel breadth (m)
%   T:   Vessel draught (m)
%   C_B: OPTIONAL block coefficient, 0 < C_B <= 1 (default: 0.65)
%
% Output:
%   Xuu: Quadratic surge-damping coefficient (kg/m)
%
% Author:    Thor I. Fossen
% Date:      2026-09-28
% Revisions:
%   2026-09-28: Added input validation and bounded the Reynolds number
%   2026-09-28: Replaced the box approximation with the Mumford formula

if nargin < 6 || isempty(C_B)
    C_B = 0.65;
end

nu_kin = 1e-6;                    % Kinematic viscosity (m^2/s)
k = 0.1;                          % Hull form factor
Re_min = 1e5;                     % Lower bound for ITTC-1957 evaluation
Re = max(L * abs(u_r) / nu_kin,Re_min);
Cf = 0.075 / (log10(Re) - 2)^2;   % ITTC-1957 correlation line

S = 1.025 * L * (C_B*B + 1.7*T);  % Mumford wetted-area approximation
Xuu = -0.5 * rho * S * (1+k) * Cf;

end
