function R = Rll(l,mu)
% R = Rll(l,mu) computes the Euler angle rotation matrix Rll 
% between ECEF and NED (l = longitude, mu = lattitude)
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 

cl  = cos(l);
sl  = sin(l);
cmu = cos(mu);
smu = sin(mu);
 
R = [...
   -cl*smu  -sl   -cl*cmu
   -sl*smu   cl   -sl*cmu
    cmu      0    -smu    ];
