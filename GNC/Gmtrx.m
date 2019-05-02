function G = Gmtrx(nabla,A_wp,GMT,GML,LCF,r_p)
% G = GMTRX(nabla,A_wp,GMT,GML,LCF,r_gp) computes the 6x6 system spring stiffness matrix G
% about an arbitrarily point P for a floating vessel (small roll and pitch angles).
% For submerged vessels, see gvect.m
% 
% Inputs:  nabla: deplasement
%          Awp: water plane area
%          GMT, GML: transverse/longitudinal metacentric heights
%          LCF = x coordinate from CO to CF, positive forwards
%                negative for conventional ships
%          r_p = [x_p y_p z_p]': location of P with respect to CO,  
%                 use  r_p = [0, 0, 0]' for CO midships
%
% Author:     Thor I. Fossen
% Date:       14th June 2001
% Revisions:  26th June 2002,  variable Awp was replaced with A_wp 
%                              one zero in G was removed
%             20th Oct 2008,   improved documentation, use r_p for
%                              arbitrarily point
%             25 Apr 2019,     added LCF as input parameter

rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

% Location of CF
r_f = [LCF, 0, 0]';

% Values in CF
G33_CF  = rho*g*A_wp;
G44_CF  = rho*g*nabla*GMT;
G55_CF  = rho*g*nabla*GML;
G_CF = diag([0 0 G33_CF G44_CF G55_CF 0]);  
G_CO = Hmtrx(r_f)' * G_CF * Hmtrx(r_f);
G = Hmtrx(r_p)' * G_CO * Hmtrx(r_p);
