function g = gvect(W,B,theta,phi,r_bG,r_bB)
% g = gvect(W,B,theta,phi,r_bB,r_bB) computes the 6x1 vector of restoring 
% forces about an arbitrarily point CO for a submerged body. Use 
% g = gRvect(W,B,R,r_bG,r_bB) if the rotation matrix R is the desired input. 
% For floating vessels, use Gmtrx.m.
% 
% Inputs: 
%  W, B: weight and buoyancy (N)
%  phi,theta: roll and pitch angles (rad)
%  r_bG = [xG yG zG]: location of the CG with respect to the CO (m)
%  r_bB = [xB yB zB]: location of the CB with respect to th CO (m)
%
% Author:    Thor I. Fossen
% Date:      14 Jun 2001
% Revisions: 20 Oct 2008 new description
%            22 Sep 2013 corrected sign error on last row (Mohammad Khani)
%            16 Dec 2021 minor updates of the documentation

sth  = sin(theta); cth  = cos(theta);
sphi = sin(phi);   cphi = cos(phi);

g = [...
   (W-B) * sth
  -(W-B) * cth * sphi
  -(W-B) * cth * cphi
  -(r_bG(2)*W-r_bB(2)*B) * cth * cphi + (r_bG(3)*W-r_bB(3)*B) * cth * sphi
   (r_bG(3)*W-r_bB(3)*B) * sth        + (r_bG(1)*W-r_bB(1)*B) * cth * cphi
  -(r_bG(1)*W-r_bB(1)*B) * cth * sphi - (r_bG(2)*W-r_bB(2)*B) * sth       ];

end
