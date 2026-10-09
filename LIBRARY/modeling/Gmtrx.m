function G = Gmtrx(nabla,A_wp,GMT,GML,x_F,r_bP)
% G = Gmtrx(nabla,A_wp,GMT,GML,x_F,r_bP) computes the 6x6 system spring 
% stiffness matrix G about an arbitrarily point P for a floating vessel 
% (small roll and pitch angles). For submerged vessels, use gvect.m or
% gRvect.m.
% 
% Inputs:  
%  nabla: Volume displacement [m^3]
%  Awp: Waterplane area [m^2]
%  GMT, GML: Transverse/longitudinal metacentric heights [m]
%  x_F = x coordinate from the CO to the CF, positive forwards negative for
%        conventional ships [m]
%  r_bP = [x_P y_P z_P]': location of the point P with respect to the CO [m],  
%          use r_bP = [0, 0, 0]' for P = CO 
%
% Author:     Thor I. Fossen
% Date:       14 Jun 2001
% Revisions:  26 Jun 2002 variable Awp was replaced by A_wp, one zero
%                         in G was removed (bugfix)
%             20 Oct 2008 updated the documentation, start using r_p for an
%                         arbitrarily point
%             25 Apr 2019 added LCF as input parameter
%             16 Dec 2021 minor updates of the documentation

constant = mssConstants();
rho = constant.rho_water;
g = constant.g;

% Location of the center of flotation (CF)
r_bF = [x_F, 0, 0]';

% Hydrostatic quantities expressed in the CF
G33_CF  = rho * g * A_wp;
G44_CF  = rho * g * nabla * GMT;
G55_CF  = rho * g * nabla * GML;
G_CF = diag([0 0 G33_CF G44_CF G55_CF 0]);  
G_CO = Hmtrx(r_bF)' * G_CF * Hmtrx(r_bF);
G = Hmtrx(r_bP)' * G_CO * Hmtrx(r_bP);

end
