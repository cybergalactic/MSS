function angle = ssa(angle,unit)
% SSA maps an angle to its smallest signed (principal) value.
%
% Examples:
%   angle = ssa(angle)       maps an angle in rad to [-pi, pi)
%   angle = ssa(angle,'deg') maps an angle in deg to [-180, 180)
%
% In feedback control systems and state estimators, angular differences
% should be mapped to [-pi, pi) or [-180, 180) to obtain the shortest
% signed angular difference and avoid artificial 2*pi or 360-deg jumps
% in feedback errors and innovations.
%
% Note that in some languages (C, C++, C#, JavaScript), the remainder
% operator may return a value with the same sign as the dividend.
% In this case, use a modulo function defined by
%
%   mod(x,y) = x - floor(x/y) * y
%
% For the Unity game engine, use Mathf.DeltaAngle.
%
% Author:     Thor I. Fossen
% Date:       2018-09-21
% Revisions:  
%   2020-03-04  Default rad, optional argument for degrees.
%   2024-12-04  Allow rad as optional argument. Throw error for invalid
%               unit. Author: Tor Børve Rasmussen.

if nargin == 1 || strcmp(unit, 'rad')
    angle = mod( angle + pi, 2 * pi ) - pi; 
elseif strcmp(unit,'deg')
    angle = mod( angle + 180, 360 ) - 180; 
else
    error("MSS:InvalidArgument", "Invalid unit argument: \'%s\'. " + ...
        "Unit must be either \'deg\' or \'rad\'", unit);
end
    
