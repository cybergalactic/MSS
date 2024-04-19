function rangeCheck(x, lower, upper)
% rangeCheck - Compatibel with MATLAB and GNU Octave (www.octave.org)
% 
% Validates that a value or values are within a specified range.
% rangeCheck(x, lower, upper) throws an error if any element of x
% is not within the range defined by lower and upper inclusive. This
% function is used to enforce that variables meet expected range conditions.
%
% Inputs:
%   x:     Numeric scalar or array to be checked.
%   lower: Numeric scalar specifying the minimum allowable value of x.
%   upper: Numeric scalar specifying the maximum allowable value of x.
%
%  Outputs:
%   None
%
%  Example:
%   rangeCheck(500, 0, 1525) % No error expected
%   rangeCheck([-10, 100, 200], 0, 150) % Error expected
%
% See also: mustBeInRange
%
% Author:     Thor I. Fossen
% Date:       2024-04-19

if any(x < lower) || any(x > upper)
    error('Value must be in the range [%d, %d].', lower, upper);
end

end
