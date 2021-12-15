function T = Tquat(u)
% Tq = Tquat(q) computes the quaternion transformation matrix Tq of 
%   dimension 4 x 3 for attitude such that q_dot = Tq * w 
% Tw = Tquat(w) computes the quaternion transformation matrix Tw of
%   dimension 4 x 4 for attitude such that q_dot = Tw * q
%
% The standard inout is the unit quaternion q = [ eta eps1 eps2 eps3 ]' of 
% dimenison 4. If the angular rates w ≈ [ p q r ]' of dimension 3 are used 
% as input, the alternative matrix Tw is returned.
%
% Author:    Thor I. Fossen
% Date:      10 Sep 2010 - core algorithm   q_dot = Tquat(q) w
% Revisions: 13 Dec 2021 - added option for q_dot = Tquat(w) q
 
if (length(u) == 4)                % T = Tq and q = u
    
    eta  = u(1); eps1 = u(2); eps2 = u(3); eps3 = u(4);
    
    T = 0.5 * [...
        -eps1 -eps2 -eps3        
         eta  -eps3  eps2
         eps3  eta  -eps1
        -eps2  eps1  eta   ];
    
elseif (length(u) == 3)            % T = Tw and w = u
    
    w = [u(1) u(2) u(3)]';
    
    T = 0.5 * [...
        0  -w'
        w  -Smtrx(w) ];
    
else
    
    error('input must be of dim. 4 (unit quaternion) or dim. 3 (angular rates)') 
    
end
