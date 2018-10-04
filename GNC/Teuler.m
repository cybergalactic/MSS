function Tzyx = Teuler(phi,theta,psi)
% Tzyx = Teuler(phi,theta,psi) computes the Euler angle
% transformation matrix for attitude
%
% Author:   Thor I. Fossen
% Date:     10th September 2010
% Revisions: 

cphi = cos(phi);
sphi = sin(phi);
cth  = cos(theta);
sth  = sin(theta);
 
J1 = Rzyx(phi,theta,psi);
 
if cth==0, error('Tzyx is singular for theta = +-90 degrees'); end
 
Tzyx = [...
      1  sphi*sth/cth  cphi*sth/cth;
      0  cphi          -sphi;
      0  sphi/cth      cphi/cth ];