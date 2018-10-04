function G = Gmtrx(nabla,A_wp,GMT,GML,r_p)
% G = GMTRX(nabla,A_wp,GMT,GML,r_gp) computes the 6x6 system spring stiffness matrix G
% about an arbitrarily point P for a floating vessel (small roll and pitch angles).
% For submerged vessels, see gvect.m
% 
% Inputs:  nabla: deplasement
%          Awp: water plane area
%          GMT, GML: transverse/longitudinal metacentric heights
%          r_p = [x_p y_p z_p]': location of P with respect to CO
%
% Author:     Thor I. Fossen
% Date:       14th June 2001
% Revisions:  26th June 2002,  variable Awp was replaced with A_wp 
%                              one zero in G was removed
%             20th Oct 2008,   improved documentation, use r_p for
%                              arbitrarily point

rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

Zz     = -rho*g*A_wp;
Kphi   = -rho*g*nabla*GMT;
Mtheta = -rho*g*nabla*GML;
G_CO   = diag([0 0 -Zz -Kphi -Mtheta 0]);  % assumes that CO = CF
G = Hmtrx(r_p)' * G_CO * Hmtrx(r_p);