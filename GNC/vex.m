function a = vex(S)
% a = vex(S) computes a = vex(S(a)) where S is a skew-symmetric matrix:
%
% S = [    0  -a(3)   a(2)
%        a(3)     0  -a(1)
%       -a(2)   a(1)   0 ];
% 
% See also: S = Smtrx(a) 
%
% Author:   Thor I. Fossen
% Date:      14th June 2001
% Revisions: 21 April 2019 - minor updates
 

if isequal(S,-S')
    a = [ S(3,2) S(1,3) S(2,1)]';
else
    error('The input argument is not a 3 x 3 skew-symmetric matrix')
end
