function [A11,ratio] = addedMassSurge(m,L,rho)
% A11 = addedMassSurge(m,L,rho) approximates the added mass in surge by
% the formula of Söding (1982):
%
%   A11 = -Xudot = 2.7 * rho * nabla^(5/3) / L^2
%
% Output:  
%   A11: added mass in surge due to a linear acceleration in surge (kg)
%   ratio: optional output for computation of the mass ratio A11/m,
%          typical values are ratios < 0.20.
%
% Inputs:   
%   m: ship/vehicle mass (kg)
%   L: ship/vehicle length (m)
%   rho: density of water, default value 1025 (kg/m3)
%
% Examples:
%   A11 = addedMassSurge(m,L)         - use default rho = 1025 (kg/m3)
%   A11 = addedMassSurge(m,L,rho)
%   [A11,ratio] = addedMassSurge(m,L,rho)
%
% Reference: H. Söding (1982). Prediction of Ship Steering Capabilities. 
%   Schiffstechnik, 3-29.
%  
% Author:    Thor I. Fossen
% Date:      17 Dec 2021
% Revisions: 

if (nargin == 2)
    rho = 1025;                         % default density of water (kg/m3)
end

nabla = m / rho;                        % volume displacement
A11 = 2.7 * rho * nabla^(5/3) / L^2;    % added mass A11
ratio = A11 / m;

end
