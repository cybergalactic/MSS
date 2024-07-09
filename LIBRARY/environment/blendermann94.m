function [tau_w, CX, CY, CK, CN] = blendermann94(...
    gamma_r, V_r, AFw, ALw, sH, sL, Loa, vessel_no)
% blendermann94 is compatible with MATLAB and GNU Octave (www.octave.org).
% The function [tau_w,CX,CY,CK,CN] = blendermann94(gamma_r,V_r,AFw,ALw,sH,
% sL,Loa,vessel_no) returns the the wind force/moment vector 
% tau_w = [tausX, tauY, tauN] and the optionally wind coeffisients CX, CY,
% CK, and CN for merchant ships using the formulas of Blendermann (1994). 
%
% INPUTS:
% gamma_r   : Relative wind angle (rad)
% V_r       : Relative wind speed (m/s)
% ALw       : Lateral projected area (m^2)
% AFw       : Frontal projected area (m^2)
% sH        : Horizontal distance to centroid of ALw (from main section)
% sL        : Vertical distance to centroid of ALw (from water line)
% Loa       : Length overall (m)
% vessel_no :
%   1. Car carrier	
%   2. Cargo vessel, loaded	
%   3. Cargo vessel, container on deck	
%   4. Container ship, loaded
%   5. Destroyer
%   6. Diving support vessel
%   7. Drilling vessel
%   8. Ferry	
%   9. Fishing vessel	
%   10. Liquified natural gas tanker	
%   11. Offshore supply vessel	
%   12. Passenger liner	
%   13. Research vessel
%   14. Speed boat	
%   15. Tanker, loaded	
%   16. Tanker, in ballast	
%   17. Tender	
%
% References:
%
%   W. Blendermann (1994). Parameter Identiﬁcation of Wind Loads on Ships. 
%     Journal of Wind Engineering and Industrial Aerodynamics 51, pp. 339–351.
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%    Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% See also: isherwood72
%
% Author:    Thor I. Fossen
% Date:      2008-11-20
% Revisions: 

if nargin~=8, error('The number of inputs must be 8');end

% Conversions and constants
rho_a = 1.224;                                  % Density of air at 20 C

% BDATA = [CD_t	CD_l_AF(0) CD_l_AF(pi) delta kappa (Fossen 2021, Ch. 10.1.2)
BDATA = [...
0.95	0.55	0.60	0.80	1.2
0.85	0.65	0.55	0.40	1.7
0.85	0.55	0.50	0.40	1.4
0.90	0.55	0.55	0.40	1.4
0.85	0.60	0.65	0.65	1.1
0.90	0.60	0.80	0.55	1.7
1.00	0.5*(0.70+1.00)	0.5*(0.75+1.10)	0.10	1.7
0.90	0.45	0.50	0.80	1.1
0.95	0.70	0.70	0.40	1.1
0.70	0.60	0.65	0.50	1.1
0.90	0.55	0.80	0.55	1.2
0.90	0.40	0.40	0.80	1.2
0.85	0.55	0.65	0.60	1.4
0.90	0.55	0.60	0.60	1.1
0.70	0.90	0.55	0.40	3.1
0.70	0.75	0.55	0.40	2.2
0.85	0.55	0.55	0.65	1.1 ];

CDt             = BDATA(vessel_no,1);
CDl_AF_bow      = BDATA(vessel_no,2);
CDl_AF_stern    = BDATA(vessel_no,3);
delta           = BDATA(vessel_no,4);
kappa           = BDATA(vessel_no,5);

Hm  = ALw/Loa; 

% Two cases for CDl
for i = 1:length(gamma_r)
    if abs(gamma_r(i))<= pi/2
        CDlAF(i,1) = CDl_AF_bow;
    else
        CDlAF(i,1) = CDl_AF_stern;
    end
end

% Wind coefficients
CDl = CDlAF * AFw / ALw;
den = 1 - 0.5 * delta * (1-CDl/CDt) .* sin(2 * gamma_r).^2;

CX = -CDlAF .* cos(gamma_r) ./ den;
CY =  CDt * sin(gamma_r) ./ den;
CK =  kappa * (sH/Hm) * CY;
CN =  ( sL/Loa - 0.18 * (gamma_r - pi/2) ) .* CY;

% Wind forces and moment
tauX = 0.5 * CX * rho_a * V_r^2 * AFw;
tauY = 0.5 * CY * rho_a * V_r^2 * ALw;
tauN = 0.5 * CN * rho_a * V_r^2 * ALw * Loa;

tau_w = [tauX,tauY,tauN]';