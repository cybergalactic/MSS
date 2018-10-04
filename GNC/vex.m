function a = vex(S)
% a = vex(S) computes the inverse of the skew-symmetric matrix:
%
% S = [    0  -a(3)   a(2)
%        a(3)     0  -a(1)
%       -a(2)   a(1)   0 ];
% 
% See also: S = Smtrx(a) 
%
% Author:   Thor I. Fossen
% Date:     6th August 2011
% Revisions: 
 
a = zeros(3,1);
if S+S'== 0
    a(1) = S(3,2);
    a(2) = S(1,3);
    a(3) = S(2,1);
else
    error('The input argument is not a 3 x 3 skew-symmetric matrix')
end
