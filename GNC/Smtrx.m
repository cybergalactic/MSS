function S = Smtrx(a)
% S = Smtrx(a) computes the 3x3 vector skew-symmetric matrix S(a) = -S(a)'.
% The corss product satisfies: a x b = S(a)b. The inverse can be computed
% as: vex(S(a)) = a 
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 6th August 2011 - updated documentation
 
S = [    0  -a(3)   a(2)
      a(3)     0   -a(1)
     -a(2)   a(1)     0 ];
 