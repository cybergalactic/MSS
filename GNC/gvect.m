function g = gvect(W,B,theta,phi,r_bg,r_bb)
% g = gvect(W,B,theta,phi,r_bg,r_bb) computes the 6x1 vector of restoring 
% forces about an arbitrarily point CO for a submerged body. For floating 
% vessels, use Gmtrx.m
% 
% Inputs:  W, B: weight and buoyancy
%          phi,theta: roll and pitch angles
%          r_bg = [x_g y_g z_g]: location of the CG with respect to the CO
%          r_bb = [x_b y_b z_b]: location of the CB with respect to th CO
%
% Author:    Thor I. Fossen
% Date:      14th June 2001
% Revisions: 20 oct 2008 improved documentation
%            22 sep 2013 corrected sign error on last row (Mohammad Khani)
%            24 Apr 2021, minor updates

sth  = sin(theta); cth  = cos(theta);
sphi = sin(phi);   cphi = cos(phi);

g = [...
   (W-B) * sth
  -(W-B) * cth * sphi
  -(W-B) * cth * cphi
  -(r_bg(2)*W-r_bb(2)*B) * cth * cphi + (r_bg(3)*W-r_bb(3)*B) * cth * sphi
   (r_bg(3)*W-r_bb(3)*B) * sth        + (r_bg(1)*W-r_bb(1)*B) * cth * cphi
  -(r_bg(1)*W-r_bb(1)*B) *cth * sphi  - (r_bg(2)*W-r_bb(2)*B) * sth       ];
