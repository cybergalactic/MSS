function [A11, ratio] = addedMassSurge(m, L, rho)
% addedMassSurge is compatible with MATLAB and GNU Octave (www.octave.org). 
% The function A11 = addedMassSurge(m,L,rho) approximates the added mass 
% in surge bythe formula of Söding (1982):
%
%   A11 = -Xudot = 2.7 * rho * nabla^(5/3) / L^2
%
% Inputs:
%   m:     mass of the ship or vehicle (kg)
%   L:     length of the ship or vehicle (m)
%   rho:   density of water (kg/m3), default value is 1025 kg/m3
%
% Outputs:
%   A11:   added mass in surge (kg)
%   ratio: optional output, ratio of added mass to actual mass (A11/m)
%
% Example usage:
%   A11 = addedMassSurge(m, L)         - use default rho = 1025 (kg/m3)
%   A11 = addedMassSurge(m, L, rho)
%   [A11,ratio] = addedMassSurge(m, L, rho)
%
% Reference: H. Söding (1982). Prediction of Ship Steering Capabilities. 
%   Schiffstechnik, 3-29.
%  
% Author:    Thor I. Fossen
% Date:      2021-12-17
% Revisions: 

if nargin == 2      % Check if the density of water is not specified
    rho = 1025;     % Default density of water (kg/m^3)
end

nabla = m / rho;                       % Volume displacement (m^3)
A11 = 2.7 * rho * nabla^(5/3) / L^2;   % Added mass in surge (kg)
ratio = A11 / m;                       % Ratio of added mass to actual mass

end
