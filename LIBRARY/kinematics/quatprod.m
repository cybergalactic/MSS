function q = quatprod(q1, q2)
% quatprod is compatible with MATLAB and GNU Octave (www.octave.org).
% This function calculates the quaternion product q = q1 o q2 of two 
% unit quaternions q1 and q2. The quaternions are assumed to be in the 
% scalar-first format, where the first element is the real part, and the 
% remaining three elements are the imaginary parts.
%
% Author: Thor I. Fossen
% Date: 2010-09-10
% Revisions:

q1 = q1(:);  % Ensure q1 is a column vector
q2 = q2(:);  % Ensure q2 is a column vector

eta1 = q1(1); 
eps1 = q1(2:4);
eta2 = q2(1); 
eps2 = q2(2:4);

q = [ eta1 * eta2 - eps1' * eps2;
      eta2 * eps1 + eta1 * eps2 + cross(eps1, eps2) ];
end