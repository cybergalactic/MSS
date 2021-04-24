function G = Gmtrx(nabla,A_wp,GMT,GML,LCF,r_bp)
% G = Gmtrx(nabla,A_wp,GMT,GML,LCF,r_bp) computes the 6x6 system spring 
% stiffness matrix G about an arbitrarily point P for a floating vessel 
% (small roll and pitch angles). For submerged vessels, use gvect.m
% 
% Inputs:  nabla: deplasement
%          Awp: waterplane area
%          GMT, GML: transverse/longitudinal metacentric heights
%          LCF = x coordinate from the CO to the CF, positive forwards
%                negative for conventional ships
%          r_bp = [x_p y_p z_p]': location of P with respect to the CO,  
%                 use  r_bp = [0, 0, 0]' for CO midships
%
% Author:     Thor I. Fossen
% Date:       14 June 2001
% Revisions:  26 June 2002,  variable Awp was replaced with A_wp 
%                            one zero in G was removed
%             20 Oct 2008,   improved documentation, use r_p for
%                            arbitrarily point
%             25 Apr 2019,   added LCF as input parameter
%             24 Apr 2021,   improved documentation

rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

% Location of the CF
r_bf = [LCF, 0, 0]';

% Hydorstatic quantities expressed in the CF
G33_CF  = rho * g * A_wp;
G44_CF  = rho * g * nabla * GMT;
G55_CF  = rho * g * nabla * GML;
G_CF = diag([0 0 G33_CF G44_CF G55_CF 0]);  
G_CO = Hmtrx(r_bf)' * G_CF * Hmtrx(r_bf);
G = Hmtrx(r_bp)' * G_CO * Hmtrx(r_bp);
