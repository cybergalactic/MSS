function g = gvect(W,B,theta,phi,r_g,r_b)
% g = GVECT(W,B,theta,phi,r_g,r_b) computes the 6x1 vector of restoring 
% forces about an arbitrarily point CO for a submerged body. For floating 
% vessels, see Gmtrx.m
% 
% Inputs:  W, B: weight and buoyancy
%          phi,theta: roll and pitch angles
%          r_g = [x_g y_g z_g]: location of CG with respect to CO
%          r_b = [x_b y_b z_b]: location of CB with respect to CO
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 20 oct 2008 improved documentation
%            22 sep 2013 corrected sign error on last row (Mohammad Khani)

sth  = sin(theta); cth  = cos(theta);
sphi = sin(phi);   cphi = cos(phi);

g = [...
   (W-B)*sth
  -(W-B)*cth*sphi
  -(W-B)*cth*cphi
  -(r_g(2)*W-r_b(2)*B)*cth*cphi + (r_g(3)*W-r_b(3)*B)*cth*sphi
   (r_g(3)*W-r_b(3)*B)*sth      + (r_g(1)*W-r_b(1)*B)*cth*cphi
  -(r_g(1)*W-r_b(1)*B)*cth*sphi - (r_g(2)*W-r_b(2)*B)*sth      ];
