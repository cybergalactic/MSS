function T = Tzyx(phi,theta)
% T = Tzyx(phi,theta) computes the Euler angle
% transformation matrix T for attitude using the zyx convention
%
% Author:   Thor I. Fossen
% Date:     4th August 2011
% Revisions: 

cphi = cos(phi);
sphi = sin(phi);
cth  = cos(theta);
sth  = sin(theta);
 
if cth==0, error('Tzyx is singular for theta = +-90 degrees'); end
 
T = [...
      1  sphi*sth/cth  cphi*sth/cth;
      0  cphi          -sphi;
      0  sphi/cth      cphi/cth ];
