function g = gRvect(W,B,R,r_bG,r_bB)
% g = gRvect(W,B,R,r_bG,r_bB) computes the 6x1 vector of restoring 
% forces about an arbitrarily point CO for a submerged body using the 
% rotation matrix R as input. Use g = gvect(W,B,theta,phi,r_bG,r_bB) for 
% Euler angle inputs. For floating vessels, use Gmtrx.m. The following 
% examples show how to use gvect.m and gRvect.m.
%
% Euler angles: 
%  R = Rzyx(phi,theta,psi)
%  g = gRvect(W,B,R,r_bG,r_bB)            input: rotation matrix R
%  g = gvect(W,B,theta,phi,r_bG,r_bB)     input: Euler angles phi, theta
% Unit quaternions:
%  Rq = Rquat(q)
%  g = gRvect(W,B,Rq,r_bG,r_bB)           input: rotation matrix Rq
%
% Inputs: 
%  W, B: weight and buoyancy (kg)
%  R: rotation matrix Rzyx (Euler angles) or Rquat (unit quaternions)
%  r_bG = [xG yG zG]: location of the CG with respect to the CO (m)
%  r_bB = [xB yB zB]: location of the CB with respect to th CO (m)
%
% Author:    Thor I. Fossen
% Date:      16 Dec 2021
% Revisions: 

g = [...
   -(W-B) * R(3,1)
   -(W-B) * R(3,2)
   -(W-B) * R(3,3)
   -(r_bG(2)*W - r_bB(2)*B) * R(3,3) + (r_bG(3)*W - r_bB(3)*B) * R(3,2)
   -(r_bG(3)*W - r_bB(3)*B) * R(3,1) + (r_bG(1)*W - r_bB(1)*B) * R(3,3)
   -(r_bG(1)*W - r_bB(1)*B) * R(3,2) + (r_bG(2)*W - r_bB(2)*B) * R(3,1) ];

end
   
