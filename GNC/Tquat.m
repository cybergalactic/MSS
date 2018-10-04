function Tq = Tquat(q)
% Tq = TQUAT(q) computes the quaternion transformation matrix for attitude
%
% Author:   Thor I. Fossen
% Date:     10th September 2010
% Revisions: 

eta  = q(1); eps1 = q(2); eps2 = q(3); eps3 = q(4); 
 
Tq = 0.5*[...
   -eps1 -eps2 -eps3        
    eta  -eps3  eps2
    eps3  eta  -eps1
   -eps2  eps1  eta   ];
