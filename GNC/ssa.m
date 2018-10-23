function angle = ssa(angle,unit)
% SSA is the "smallest signed angle" or the smallest difference between two
% angles. Examples:
%  
% angle = ssa(angle,'rad') maps an angle in rad to the interval [-pi pi) 
% angle = ssa(angle,'deg') maps an angle in deg to the interval [-180, 180)
%
% For feedback control systems and state estimators used to control the 
% attitude of vehicles, the difference of two angles should always be
% mapped to [-pi pi) or [-180, 180] to avoid step inputs/discontinuties.           
%
% Notice that in many languages, the modulus operator mod(x,y) returns a
% value with the same sign as x:
%
% - C, C++, C# and JavaScript: Use a custom mod function e.g.: 
%   mod(x, y) = x - floor(x/y) * y;
% - Unity game engine: Use Mathf.DeltaAngle
%
% Author:     Thor I. Fossen
% Date:       2018-09-21
% Revisions:  
%__________________________________________________________________________

if strcmp(unit,'rad')
    angle = mod( angle + pi, 2 * pi ) - pi; 
elseif strcmp(unit,'deg')
    angle = mod( angle + 180, 360 ) - 180; 
else
    error('- Wrong unit in ssa(angle,unit). Use ''rad'' or ''deg''.');
end
    
