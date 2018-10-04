function [l,mu,h] = ecef2llh(x,y,z)
% [l,mu,h] = ECEF2LLH(x,y,z) computes the longitude l (rad),
% latitude mu (rad) and height h (m) from the ECEF positions (x,y,z)
%
% Author:   Thor I. Fossen
% Date:     7th June 2001
% Revisions: 1st September 2002, atan2(y/x) replaced by atan2(y,x) 
%            2nd September 2002, new output argument for height h was added
%            27th January, 2003, angle outputs are defined in rad

r_e = 6378137;          % WGS-84 data
r_p = 6356752;
e = 0.08181979099211;
l = atan2(y,x);
eps = 1;
tol = 1e-10;
p = sqrt(x^2+y^2);
mu = atan(z/(p*(1-e^2)));

while (eps > tol),
   N = r_e^2/sqrt(r_e^2*cos(mu)^2+r_p^2*sin(mu)^2);
   h = p/cos(mu)-N;
   mu0 = mu;
   mu = atan(z/(p*(1-e^2*N/(N+h))));
   eps = abs(mu-mu0);
end
